/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2022 OpenCFD Ltd.
    Modified code Copyright (C) 2019-2020 DLR
    Modified code Copyright (C) 2018, 2021 Johan Roenby
    Modified code Copyright (C) 2024 Dezhi Dai, ANL
-------------------------------------------------------------------------------
License
    This file is part of geometricVofExt, which is a geometric VOF extension
    to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solveVofEqu.H"
#include "primitiveMeshTools.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "volPointInterpolation.H"
#include "fvcSurfaceIntegrate.H"
#include "cellCellStencilObject.H"
#include "fvcGrad.H"
#include "cellSet.H"
#include "meshTools.H"
#include "OFstream.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::geometricVofExt::SimPLIC::advection::typeName = "advection";


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::geometricVofExt::SimPLIC::advection::extendMarkedCells(bitSet& markedCell) const
{
    // Mark faces using any marked cell
    bitSet markedFace(mesh_.nFaces());

    for (const label celli : markedCell)
    {
        markedFace.set(mesh_.cells()[celli]);  // set multiple faces
    }

    syncTools::syncFaceList(mesh_, markedFace, orEqOp<unsigned int>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
            markedCell.set(mesh_.faceNeighbour()[facei]);
        }
    }
    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); ++facei)
    {
        if (markedFace.test(facei))
        {
            markedCell.set(mesh_.faceOwner()[facei]);
        }
    }
}


void Foam::geometricVofExt::SimPLIC::advection::timeIntegratedFlux()
{
    // Get time step
    const scalar dt(mesh_.time().deltaTValue());

    // Create object for interpolating velocity to interface centres
    interpolationCellPoint<vector> UInterp(U_);
    //interpolationCellPointFace<vector> UInterp(U_);

    // Clear out the data for re-use
    clearBoundaryData();

    // Get necessary references
    const scalarField& phiIn(phi_.primitiveField());
    const scalarField& magSfIn(mesh_.magSf().primitiveField());
    scalarField& dVfIn(dVf_.primitiveFieldRef());

    // Get necessary mesh data
    const cellList& cellFaces(mesh_.cells());
    const labelList& own(mesh_.faceOwner());

    // Loop through all mixed cells
    const DynamicList<label>& mixedCells(reconstructor_.mixedCells());
    const DynamicList<label>& cellStatus(reconstructor_.cellStatus());
    const volVectorField& interfaceN(reconstructor_.interfaceN());
    const volScalarField& interfaceD(reconstructor_.interfaceD());
    const volVectorField& interfaceC(reconstructor_.interfaceC());
    forAll(mixedCells, i)
    {
        if (cellStatus[i] != 0)
        {
            continue;
        }

        const label cellI(mixedCells[i]);
        const vector& normalI(interfaceN.primitiveField()[cellI]);
        const scalar distanceI(interfaceD.primitiveField()[cellI]);
        const point& centreI(interfaceC.primitiveField()[cellI]);

        // Get the speed of the plicInterface by interpolating velocity and
        // dotting its normal vector
        const scalar Un0(UInterp.interpolate(centreI, cellI) & normalI);

        // Estimate time integrated flux through each downwind face
        // Note: looping over all cell faces - in reduced-D, some of
        //       these faces will be on empty patches
        const cell& cellIFaces(cellFaces[cellI]);
        for (const label facei : cellIFaces)
        {
            if (mesh_.isInternalFace(facei))
            {
                bool isDownwindFace(false);

                if (cellI == own[facei])
                {
                    if (phiIn[facei] >= 0.0)
                    {
                        isDownwindFace = true;
                    }
                }
                else
                {
                    if (phiIn[facei] < 0.0)
                    {
                        isDownwindFace = true;
                    }
                }

                if (isDownwindFace)
                {
                    dVfIn[facei] = cutFace_.timeIntegratedFaceFlux
                    (
                        facei,
                        normalI,
                        distanceI,
                        Un0,
                        dt,
                        phiIn[facei],
                        magSfIn[facei]
                    );
                }
            }
            else
            {
                bsFaces_.append(facei);
                bsn0_.append(normalI);
                bsD0_.append(distanceI);
                bsUn0_.append(Un0);
            }
        }
    }

    // Get references to boundary fields
    const polyBoundaryMesh& boundaryMesh(mesh_.boundaryMesh());
    const surfaceScalarField::Boundary& phib(phi_.boundaryField());
    const surfaceScalarField::Boundary& magSfb(mesh_.magSf().boundaryField());
    surfaceScalarField::Boundary& dVfb(dVf_.boundaryFieldRef());
    const label nInternalFaces(mesh_.nInternalFaces());

    // Loop through boundary surface faces
    forAll(bsFaces_, i)
    {
        // Get boundary face index (in the global list)
        const label faceI(bsFaces_[i]);
        const label patchI(boundaryMesh.patchID()[faceI - nInternalFaces]);
        const label start(boundaryMesh[patchI].start());

        if (!phib[patchI].empty())
        {
            const label patchFaceI(faceI - start);
            const scalar phiP(phib[patchI][patchFaceI]);

            if (phiP >= 0.0)
            {
                const scalar magSf(magSfb[patchI][patchFaceI]);

                dVfb[patchI][patchFaceI] = cutFace_.timeIntegratedFaceFlux
                (
                    faceI,
                    bsn0_[i],
                    bsD0_[i],
                    bsUn0_[i],
                    dt,
                    phiP,
                    magSf
                );

                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(faceI);
            }
        }
    }

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_);
}


void Foam::geometricVofExt::SimPLIC::advection::setDownwindFaces
(
    const label cellI,
    DynamicList<label>& downwindFaces
) const
{
    // Get necessary mesh data and cell information
    const labelList& own(mesh_.faceOwner());
    const cellList& cells(mesh_.cells());
    const cell& c(cells[cellI]);

    downwindFaces.clear();

    // Check all faces of the cell
    for (const label facei: c)
    {
        const scalar phi(faceValue(phi_, facei));

        if (own[facei] == cellI)
        {
            if (phi >= 0)
            {
                downwindFaces.append(facei);
            }
        }
        else if (phi < 0)
        {
            downwindFaces.append(facei);
        }
    }

    downwindFaces.shrink();
}


Foam::scalar Foam::geometricVofExt::SimPLIC::advection::netFlux
(
    const label cellI,
    const surfaceScalarField& dVf
) const
{
    scalar dV(0.0);

    // Get face indices
    const cell& c(mesh_.cells()[cellI]);

    // Get mesh data
    const labelList& own(mesh_.faceOwner());

    for (const label facei : c)
    {
        const scalar dVff(faceValue(dVf, facei));

        if (own[facei] == cellI)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


void Foam::geometricVofExt::SimPLIC::advection::applyBruteForceBounding()
{
    if (snapAlphaTol_ > 0.0)
    {
        alpha1_ = alpha1_
                * pos0(alpha1_ - snapAlphaTol_)
                * neg0(alpha1_ - (1.0 - snapAlphaTol_))
                + pos0(alpha1_ - (1.0 - snapAlphaTol_));

        alpha1_.correctBoundaryConditions();
    }

    if (clip_)
    {
        alpha1_ = min(scalar(1), max(scalar(0), alpha1_));
        alpha1_.correctBoundaryConditions();
    }
}


Foam::DynamicList<Foam::label> Foam::geometricVofExt::SimPLIC::advection::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi,
    bool returnSyncedFaces
)
{
    DynamicList<label> syncedFaces(0);
    const polyBoundaryMesh& patches(mesh_.boundaryMesh());

    if (Pstream::parRun())
    {
        DynamicList<label> neighProcs;
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (const label patchi : procPatchLabels_)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);
            const label nbrProci(procPatch.neighbProcNo());

            neighProcs.append(nbrProci);
            UOPstream toNbr(nbrProci, pBufs);

            const scalarField& pFlux(dVf.boundaryField()[patchi]);
            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            toNbr << surfCellFacesOnProcPatch << dVfPatch;
        }

        // Limited to involved neighbour procs
        pBufs.finishedNeighbourSends(neighProcs);
        //pBufs.finishedSends();

        // Receive and combine
        for (const label patchi : procPatchLabels_)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);
            const label nbrProci(procPatch.neighbProcNo());

            UIPstream fromNeighb(nbrProci, pBufs);
            List<label> faceIDs;
            List<scalar> nbrdVfs;

            fromNeighb >> faceIDs >> nbrdVfs;
            if (returnSyncedFaces)
            {
                List<label> syncedFaceI(faceIDs);
                for (label& faceI : syncedFaceI)
                {
                    faceI += procPatch.start();
                }
                syncedFaces.append(syncedFaceI);
            }

            // Combine fluxes
            scalarField& localFlux(dVf.boundaryFieldRef()[patchi]);

            forAll(faceIDs, i)
            {
                const label facei(faceIDs[i]);
                localFlux[facei] = - nbrdVfs[i];
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }

    return syncedFaces;
}


void Foam::geometricVofExt::SimPLIC::advection::checkIfOnProcPatch(const label faceI)
{
    if (!mesh_.isInternalFace(faceI))
    {
        const polyBoundaryMesh& pbm(mesh_.boundaryMesh());
        const label patchi(pbm.patchID()[faceI - mesh_.nInternalFaces()]);

        if (isA<processorPolyPatch>(pbm[patchi]) && !pbm[patchi].empty())
        {
            const label patchFacei(pbm[patchi].whichFace(faceI));
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}


void Foam::geometricVofExt::SimPLIC::advection::setProcessorPatches()
{
    // Get boundary mesh and resize the list for parallel comms
    const polyBoundaryMesh& patches(mesh_.boundaryMesh());
    surfaceCellFacesOnProcPatches_.clear();
    surfaceCellFacesOnProcPatches_.resize(patches.size());

    // Append all processor patch labels to the list
    procPatchLabels_.clear();
    forAll(patches, patchi)
    {
        if
        (
            isA<processorPolyPatch>(patches[patchi])
         && !patches[patchi].empty()
        )
        {
            procPatchLabels_.append(patchi);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::geometricVofExt::SimPLIC::advection::advection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const reconstruction& reconstructor,
    const dictionary& dict
)
:
    // External data references
    mesh_(dynamic_cast<const dynamicFvMesh&>(alpha1.mesh())),
    alpha1_(alpha1),
    alpha1In_(alpha1.ref()),
    phi_(phi),
    U_(U),
    reconstructor_(reconstructor),

    // Switches/tolerances
    nAlphaBounds_(dict.getOrDefault<label>("nAlphaBounds", 10)),
    snapAlphaTol_(dict.getOrDefault<scalar>("snapTol", 0.0)),
    clip_(dict.lookupOrDefault<bool>("clip", true)),

    // Face cutting
    cutFace_(mesh_, reconstructor_.faceFlatness()),

    // Surface liquid volume and flux
    dVf_
    (
        IOobject
        (
            "dVf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol, 0.0)
    ),
    alphaPhi_
    (
        IOobject
        (
            "alphaPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVol/dimTime, 0.0)
    ),

    // Execution performance
    advectionTime_(0.0),

    // Boundary storages
    bsFaces_(label(0.2 * mesh_.nBoundaryFaces())),
    bsn0_(bsFaces_.size()),
    bsD0_(bsFaces_.size()),
    bsUn0_(bsFaces_.size()),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0)
{
    // Prepare lists used in parallel runs
    if (Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        setProcessorPatches();
    }

    clearBoundaryData();
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //


// ************************************************************************* //