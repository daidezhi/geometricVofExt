/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2024 Dezhi Dai, Argonne National Laboratory (ANL)
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

#include "reconstruction.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "volPointInterpolation.H"
#include "cellCellStencilObject.H"
#include "fvcGrad.H"
#include "leastSquareGrad.H"
#include "zoneDistribute.H"
#include "addToRunTimeSelectionTable.H"
#include "pointLinear.H"
#include "reconstructedDistanceFunction.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{
    defineTypeNameAndDebug(reconstruction, 0);
}
}
}


const std::map<Foam::word, Foam::geometricVofExt::SimPLIC::reconstruction::validOrientationMethod_>
Foam::geometricVofExt::SimPLIC::reconstruction::validOrientationMethods_ =
{
    {"alphaGrad",          validOrientationMethod_::ALPHAGRAD},
    {"isoAlphaGrad",       validOrientationMethod_::ISOALPHAGRAD},
    {"isoRDF",             validOrientationMethod_::ISORDF},
    {"LS",                 validOrientationMethod_::ISOALPHAGRAD},
    {"RDF",                validOrientationMethod_::ISORDF}
};


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::geometricVofExt::SimPLIC::reconstruction::alreadyReconstructed
(
    bool forceUpdate
) const
{
    const Time& runTime = alpha1_.mesh().time();

    label& curTimeIndex = timeIndexAndIter_.first();
    label& curIter = timeIndexAndIter_.second();

    // Reset timeIndex and curIter
    if (curTimeIndex < runTime.timeIndex())
    {
        curTimeIndex = runTime.timeIndex();
        curIter = 0;
        return false;
    }

    if (forceUpdate)
    {
        curIter = 0;
        return false;
    }

    // Always reconstruct when subcycling
    if (runTime.subCycling() != 0)
    {
        return false;
    }

    ++curIter;
    if (curIter > 1)
    {
        return true;
    }

    return false;
}


void Foam::geometricVofExt::SimPLIC::reconstruction::calcInterfaceNFromRegAlphaGrad()
{
    vectorField& interfaceNIn(interfaceN_.primitiveFieldRef());

    interfaceN_ = -fvc::grad(alpha1_, "grad(alpha1)");
    interfaceNIn /= (mag(interfaceNIn) + SMALL);

    interfaceN_.correctBoundaryConditions();
}


void Foam::geometricVofExt::SimPLIC::reconstruction::calcInterfaceNFromIsoAlphaGrad()
{
    vectorField& interfaceNIn(interfaceN_.primitiveFieldRef());

    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());

    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);
    exchangeFields.setUpCommforZone(interfaceCell_, true);

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, mesh_.C())
    );
    Map<scalar> mapAlpha
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, alpha1_)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> alphaValues(100);

    const labelListList& stencil(exchangeFields.getStencil());

    forAll(mixedCells_, i)
    {
        const label celli(mixedCells_[i]);

        cellCentre.clear();
        alphaValues.clear();

        for (const label gblIdx : stencil[celli])
        {
            cellCentre.append
            (
                exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
            );
            alphaValues.append
            (
                exchangeFields.getValue(alpha1_, mapAlpha, gblIdx)
            );
        }

        cellCentre -= mesh_.C()[celli];

        interfaceNIn[celli] = -lsGrad.grad(cellCentre, alphaValues);
    }

    interfaceNIn /= (mag(interfaceNIn) + SMALL);

    interfaceN_.correctBoundaryConditions();
}


void Foam::geometricVofExt::SimPLIC::reconstruction::isoInterfaceGrad
(
    const volScalarField& phi,
    DynamicField<vector>& interfaceNormal
)
{
    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());
    zoneDistribute& exchangeFields = zoneDistribute::New(mesh_);

    exchangeFields.setUpCommforZone(interfaceCell_, false);

    Map<vector> mapCC
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, mesh_.C())
    );
    Map<scalar> mapAlpha
    (
        exchangeFields.getDatafromOtherProc(interfaceCell_, phi)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> phiValues(100);

    const labelListList& stencil(exchangeFields.getStencil());

    forAll(mixedCells_, i)
    {
        const label celli(mixedCells_[i]);

        cellCentre.clear();
        phiValues.clear();

        for (const label gblIdx : stencil[celli])
        {
            cellCentre.append
            (
                exchangeFields.getValue(mesh_.C(), mapCC, gblIdx)
            );
            phiValues.append
            (
                exchangeFields.getValue(phi, mapAlpha, gblIdx)
            );
        }

        cellCentre -= mesh_.C()[celli];

        interfaceNormal[i] = lsGrad.grad(cellCentre, phiValues);
    }
}


void Foam::geometricVofExt::SimPLIC::reconstruction::calcInterfaceNFromIsoRDF()
{
    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    // Initialize orientation vectors *****************************************
    vectorField& interfaceNIn(interfaceN_.primitiveFieldRef());

    Foam::reconstructedDistanceFunction RDF(mesh_);

    DynamicField<vector> interfaceNormal(0.2*mesh_.nCells());
    interfaceNormal.setSize(mixedCells_.size());

    const volVectorField& U(alpha1_.db().lookupObject<volVectorField>("U"));

    zoneDistribute& exchangeFields(zoneDistribute::New(mesh_));

    RDF.markCellsNearSurf(interfaceCell_, 1);
    const boolList& nextToInterface(RDF.nextToInterface());
    exchangeFields.updateStencil(nextToInterface);

    isoInterfaceGrad(alpha1_, interfaceNormal);
    // ************************************************************************

    // RDF iterations *********************************************************
    interfaceC_ = dimensionedVector(dimLength, Zero);
    interfaceN_ = dimensionedVector(dimless/dimLength, Zero);

    bitSet tooCoarse(mesh_.nCells(), false);

    for (label iter=0; iter < isoRDFIterations_; ++iter)
    {
        forAll(mixedCells_, i)
        {
            const label celli(mixedCells_[i]);

            interfaceN_.primitiveFieldRef()[celli] = -interfaceNormal[i].normalise(SMALL);

            if (tooCoarse.test(celli))
            {
                continue;
            }

            cellStatus_[i] = cutCell_.findSignedDistance(celli, splitWarpedFace_);
        }

        interfaceC_.correctBoundaryConditions();
        interfaceN_.correctBoundaryConditions();

        List<label> normalResCelli(mixedCells_.size());
        List<scalar> normalResNormalResidual(mixedCells_.size());
        List<scalar> normalResAvgAngle(mixedCells_.size());

        surfaceVectorField::Boundary nHatb(mesh_.Sf().boundaryField());
        nHatb *= 1.0 / (mesh_.magSf().boundaryField());

        {
            RDF.constructRDF
            (
                nextToInterface,
                interfaceC_,
                -interfaceN_,
                exchangeFields,
                false
            );

            RDF.updateContactAngle(alpha1_, U, nHatb);

            isoInterfaceGrad(RDF, interfaceNormal);

            // Calculate residual
            Map<vector> mapNormal
            (
                exchangeFields.getDatafromOtherProc(interfaceCell_, interfaceN_)
            );

            const labelListList& stencil(exchangeFields.getStencil());

            forAll(mixedCells_, i)
            {
                const label celli(mixedCells_[i]);

                if
                (
                    mag(interfaceN_[celli]) < TSMALL or
                    mag(interfaceNormal[i].normalise(SMALL)) < TSMALL
                )
                {
                    normalResCelli[i] = celli;
                    normalResNormalResidual[i] = 0.0;
                    normalResAvgAngle[i] = 0.0;
                }

                scalar avgDiffNormal(0);
                scalar maxDiffNormal(GREAT);
                scalar weight(0);

                const vector cellNormal(interfaceN_[celli]);

                forAll(stencil[celli], j)
                {
                    const label gblIdx(stencil[celli][j]);
                    vector normal
                    (
                        exchangeFields.getValue(interfaceN_, mapNormal, gblIdx)
                    );

                    if (mag(normal) >= TSMALL && j != 0)
                    {
                        vector n(normal / mag(normal));
                        scalar cosAngle(max(min((cellNormal & n), 1.0), -1.0));
                        avgDiffNormal += acos(cosAngle) * mag(normal);
                        weight += mag(normal);
                        if (cosAngle < maxDiffNormal)
                        {
                            maxDiffNormal = cosAngle;
                        }
                    }
                }

                if (weight != 0)
                {
                    avgDiffNormal /= weight;
                }
                else
                {
                    avgDiffNormal = 0;
                }

                vector newCellNormal(-interfaceNormal[i].normalise(SMALL));

                scalar normalRes((1.0 - (cellNormal & newCellNormal)));

                normalResCelli[i] = celli;
                normalResNormalResidual[i] = normalRes;
                normalResAvgAngle[i] = avgDiffNormal;
            }
        }

        label resCounter(0);
        scalar avgRes(0);
        scalar avgNormRes(0);

        forAll(normalResCelli, i)
        {
            const label celli(normalResCelli[i]);
            const scalar normalRes(normalResNormalResidual[i]);
            const scalar avgA(normalResAvgAngle[i]);

            if (avgA > 0.26 and iter > 0) // 15 deg
            {
                tooCoarse.set(celli);
            }
            else
            {
                avgRes += normalRes;
                scalar normRes(0);
                scalar discreteError(0.01 * sqr(avgA));
                if (discreteError != 0)
                {
                    normRes = normalRes / max(discreteError, isoRDFTol_);
                }
                else
                {
                    normRes = normalRes / isoRDFTol_;
                }
                avgNormRes += normRes;
                resCounter++;
            }
        }

        reduce(avgRes,sumOp<scalar>());
        reduce(avgNormRes,sumOp<scalar>());
        reduce(resCounter,sumOp<label>());

        if (resCounter == 0) // avoid division  by zero and leave loop
        {
            resCounter = 1;
            avgRes = 0;
            avgNormRes = 0;
        }

        if
        (
            (
                (
                    avgNormRes / resCounter < isoRDFRelTol_
                    or
                    avgRes/resCounter < isoRDFTol_
                )
                and
                (iter >= 1)
            )
            or
            (iter + 1 == isoRDFIterations_)
        )
        {
            break;
        }
    }
    // ************************************************************************

    interfaceN_.correctBoundaryConditions();
}


void Foam::geometricVofExt::SimPLIC::reconstruction::updateFaceFlatness()
{
    const faceList& faces(mesh_.faces());
    const pointField& points(mesh_.points());
    const vectorField& fCtrs(mesh_.faceCentres());
    const vectorField& faceAreas(mesh_.faceAreas());
    const scalarField magAreas(mag(faceAreas));

    forAll(faces, facei)
    {
        const face& fa(faces[facei]);

        if (fa.size() > 3 && magAreas[facei] > ROOTVSMALL)
        {
            const point& fc(fCtrs[facei]);
            scalar sumA(0.0);

            forAll(fa, pointi)
            {
                const point& thisPoint(points[fa[pointi]]);
                const point& nextPoint(points[fa.nextLabel(pointi)]);

                vector n = 0.5 * ((nextPoint - thisPoint) ^ (fc - thisPoint));
                sumA += mag(n);
            }

            faceFlatness_[facei] = magAreas[facei]/(sumA + ROOTVSMALL);
        }
        else
        {
            faceFlatness_[facei] = 1.0;
        }
    }

    Info<< nl
        << "SimPLIC::Mesh face flatness: min/max/avg = "
        << gMin(faceFlatness_) << "/"
        << gMax(faceFlatness_) << "/"
        << gSum(faceFlatness_ * magAreas) / gSum(magAreas) << nl
        << endl;

    if (writePlicFields_)
    {
        const polyBoundaryMesh& boundaryMesh(mesh_.boundaryMesh());
        const label nInternalFaces(mesh_.nInternalFaces());

        scalarField& zetaIn(zeta_.primitiveFieldRef());
        surfaceScalarField::Boundary& zetab(zeta_.boundaryFieldRef());

        forAll(faces, facei)
        {
            if (mesh_.isInternalFace(facei))
            {
                zetaIn[facei] = faceFlatness_[facei];
            }
            else
            {
                const label patchI(boundaryMesh.patchID()[facei - nInternalFaces]);
                const label start(boundaryMesh[patchI].start());
                const label patchFaceI(facei - start);

                zetab[patchI][patchFaceI] = faceFlatness_[facei];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::geometricVofExt::SimPLIC::reconstruction::reconstruction
(
    volScalarField& alpha1,
    const dictionary& dict
)
:
    // For sampledPlicSurface
    IOdictionary
    (
        IOobject
        (
            "reconstruction",
            alpha1.time().constant(),
            alpha1.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    // External data references
    mesh_(dynamic_cast<const fvMesh&>(alpha1.mesh())),
    alpha1_(alpha1),
    //alpha1In_(alpha1.internalField()),
    alpha1In_(alpha1.ref()),
    mixedCellTol_(dict.getOrDefault<scalar>("mixedCellTol", 1e-8)),
    writePlicFields_(dict.getOrDefault<bool>("writePlicFields", false)),
    splitWarpedFace_(dict.getOrDefault<bool>("splitWarpedFace", false)),
    orientationMethod_
    (
        dict.getOrDefault<word>("orientationMethod", word("isoAlphaGrad"))
    ),
    orientationMethodMapIter_
    (
        validOrientationMethods_.find(orientationMethod_)
    ),
    isoRDFTol_(dict.getOrDefault<scalar>("tol", 1e-6)),
    isoRDFRelTol_(dict.getOrDefault<scalar>("relTol", 0.1)),
    isoRDFIterations_(dict.getOrDefault<label>("iterations", 5)),
    mapAlphaField_(dict.getOrDefault<bool>("mapAlphaField", false)),
    faceFlatness_(mesh_.nFaces(), scalar(1)),

    // Interface fields
    interfaceN_
    (
        IOobject
        (
            "interfaceN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            writePlicFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    interfaceD_
    (
        IOobject
        (
            "interfaceD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            writePlicFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, 0.0)
    ),
    interfaceC_
    (
        IOobject
        (
            "interfaceC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            writePlicFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimLength, point::zero)
    ),
    interfaceS_
    (
        IOobject
        (
            "interfaceS",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            writePlicFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimArea, vector::zero)
    ),
    zeta_
    (
        IOobject
        (
            "zeta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            writePlicFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    ),
    cellVolume0_
    (
        IOobject
        (
            "cellVolume0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimVolume, Zero)
    ),

    interfaceCell_(alpha1_.mesh().nCells(), false),

    //
    //orientation_(orientationSchemes::New(alpha1_, interfaceN_, mixedCells_, dict)),

    // Execution performance
    reconstructionTime_(0.0),
    alphaMappingTime_(0.0),

     // Cell and face cutting
    mixedCells_(label(0.2 * mesh_.nCells())),
    cellStatus_(label(0.2 * mesh_.nCells())),
    cutFace_(mesh_, faceFlatness_),
    cutCell_
    (
        mesh_,
        faceFlatness_,
        alpha1,
        interfaceN_,
        interfaceD_,
        interfaceC_,
        interfaceS_
    ),
    timeIndexAndIter_(0, 0)
{
    clearInterfaceData();

    cellVolume0_.primitiveFieldRef() = mesh_.cellVolumes();

    if (orientationMethodMapIter_ == validOrientationMethods_.end())
    {
        FatalErrorInFunction
            << "Orientation vector calculation method '"
            << orientationMethod_
            << "' is not valid. Valid methods are" << nl
            << "(" << nl
            << "    alphaGrad" << nl
            << "    isoAlphaGrad" << nl
            << "    isoRDF " << nl
            << "    smoothAlphaGrad" << nl
            << "    smoothIsoAlphaGrad" << nl
            << "    smoothIsoRDF" << nl
            << ")"
            << abort(FatalError);
    }

    // Update face flatness field
    updateFaceFlatness();
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

void Foam::geometricVofExt::SimPLIC::reconstruction::initialize()
{
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    interfaceCell_ = false;

    clearInterfaceData();

    interfaceN_ == dimensionedVector(dimless/dimLength, Zero);
    interfaceD_ == dimensionedScalar(dimLength, Zero);
    interfaceC_ == dimensionedVector(dimLength, Zero);
    interfaceS_ == dimensionedVector(dimArea, Zero);

    if (mesh_.changing())
    {
        faceFlatness_.resize(mesh_.nFaces(), scalar(1));
        updateFaceFlatness();
    }

    if (isA<dynamicOversetFvMesh>(mesh_))
    {
        const cellCellStencilObject& overlap(Stencil::New(mesh_));
        const labelList& cellTypes(overlap.cellTypes());

        forAll(alpha1In_, celli)
        {
            if (isAMixedCell(celli) && cellTypes[celli]==cellCellStencil::CALCULATED)
            {
                mixedCells_.append(celli);
                interfaceCell_[celli] = true;
                cellStatus_.append(-100);
            }
        }
    }
    else
    {
        forAll(alpha1In_, celli)
        {
            if (isAMixedCell(celli))
            {
                mixedCells_.append(celli);
                interfaceCell_[celli] = true;
                cellStatus_.append(-100);
            }
        }
    }

    Info<< "SimPLIC::reconstruction: Number of mixed cells = "
        << returnReduce(mixedCells_.size(), sumOp<label>()) << endl;
}


void Foam::geometricVofExt::SimPLIC::reconstruction::reconstruct(bool forceUpdate)
{
    // Remember the volumes of the mesh the interface is reconstructed on.
    // The field is mapped by mesh.update(), which lets mapAlphaField()
    // identify cells created by refinement.
    cellVolume0_.primitiveFieldRef() = mesh_.cellVolumes();

    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    const scalar startTime(mesh_.time().elapsedCpuTime());

    initialize();

    if (returnReduce(mixedCells_.size(), sumOp<label>()) == 0)
    {
        return;
    }

    // Step 1: Calculate interface orientations * * * * * * * * * * * * * * * *
    switch (orientationMethodMapIter_->second)
    {
        case validOrientationMethod_::ALPHAGRAD:
            calcInterfaceNFromRegAlphaGrad();
            break;

        case validOrientationMethod_::ISOALPHAGRAD:
            calcInterfaceNFromIsoAlphaGrad();
            break;

        case validOrientationMethod_::ISORDF:
            calcInterfaceNFromIsoRDF();
            break;
    }
    //  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


    // Step 2: Interface plane locating * * * * * * * * * * * * * * * * * * * *
    forAll(mixedCells_, i)
    {
        cellStatus_[i] =
            cutCell_.findSignedDistance(mixedCells_[i], splitWarpedFace_);

        if (cellStatus_[i] == 0)
        {
            if(mag(interfaceN_.ref()[mixedCells_[i]]) <= SMALL)
            {
                interfaceCell_[mixedCells_[i]] = false;
            }
        }
        else
        {
            interfaceCell_[mixedCells_[i]] = false;

            interfaceN_.primitiveFieldRef()[mixedCells_[i]] = Zero;
            interfaceD_.primitiveFieldRef()[mixedCells_[i]] = Zero;
            interfaceC_.primitiveFieldRef()[mixedCells_[i]] = Zero;
            interfaceS_.primitiveFieldRef()[mixedCells_[i]] = Zero;
        }
    }

    interfaceD_.correctBoundaryConditions();
    interfaceC_.correctBoundaryConditions();
    interfaceS_.correctBoundaryConditions();
    //  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    reconstructionTime_ += (mesh_.time().elapsedCpuTime() - startTime);
}


void Foam::geometricVofExt::SimPLIC::reconstruction::mapAlphaField()
{
    // Only meaningful after a topology change of an AMR mesh.
    // Note: the previous check
    //     if (!isA<dynamicRefineFvMesh>(mesh_) or !isA<dynamicRefine2DFvMesh>(mesh_))
    // was always true (a mesh cannot be both types), which silently disabled
    // the geometric mapping for every mesh type.
    if (!mapAlphaField_ || !isRefineMesh() || !mesh_.topoChanging())
    {
        return;
    }

    const scalar startTime(mesh_.time().elapsedCpuTime());

    // A child cell created by refinement has V/V0 = 1/4 (2D) or 1/8 (3D),
    // where V0 is the (mapped) volume of its parent. Unchanged cells have
    // V/V0 = 1 and merged (unrefined) cells have V/V0 > 1.
    const scalar refinedVolumeRatio(0.6);

    // Reconstructed normals are unit vectors; cells without a valid interface
    // (pure cells, or mixed cells whose reconstruction failed) carry n = 0.
    const scalar validNormalMag(0.5);

    const scalarField& V(mesh_.cellVolumes());
    const scalarField& V0(cellVolume0_.primitiveField());
    const vectorField& normals(interfaceN_.primitiveField());
    const scalarField& distances(interfaceD_.primitiveField());

    label nRefinedChildren(0);
    label nMappedCells(0);
    label nMixedWithoutNormal(0);

    // Phase volume before geometric mapping (after the mesh mapper)
    const scalar phaseVolume0(gSum(alpha1In_*V));

    forAll(alpha1In_, celli)
    {
        // Only children of refined parents are mapped geometrically. Their
        // interfaceN/interfaceD were copied from the parent, so cutting all
        // children with the parent plane conserves the parent volume.
        // Unchanged and merged cells keep the (conservative) values given
        // by the mesh mapper.
        if (V[celli] > refinedVolumeRatio*V0[celli])
        {
            continue;
        }

        ++nRefinedChildren;

        // Parent had no valid interface plane: the copied alpha is kept.
        // (Previously such cells were cut with n = 0 and set to alpha = 0,
        // i.e. liquid was deleted.)
        if (mag(normals[celli]) < validNormalMag)
        {
            // Copied alpha is mixed, i.e. the parent was an interface cell,
            // but no plane is available: should not happen in normal runs
            if (isAMixedCell(celli))
            {
                ++nMixedWithoutNormal;
            }

            continue;
        }

        cutCell_.calcSubCell(celli, normals[celli], distances[celli], false);

        alpha1In_[celli] = cutCell_.volumeOfFluid();

        ++nMappedCells;
    }

    const scalar phaseVolume1(gSum(alpha1In_*V));

    alpha1_.correctBoundaryConditions();
    alpha1_.oldTime() = alpha1_;
    alpha1_.oldTime().correctBoundaryConditions();

    Info<< "SimPLIC::reconstruction: Refined child cells = "
        << returnReduce(nRefinedChildren, sumOp<label>())
        << ", geometrically mapped = "
        << returnReduce(nMappedCells, sumOp<label>())
        << ", mixed without normal = "
        << returnReduce(nMixedWithoutNormal, sumOp<label>())
        << ", phase volume change = " << (phaseVolume1 - phaseVolume0)
        << endl;

    alphaMappingTime_ += (mesh_.time().elapsedCpuTime() - startTime);
}


bool Foam::geometricVofExt::SimPLIC::reconstruction::isRefineMesh() const
{
    return
    (
        isA<dynamicRefineFvMesh>(mesh_) || isA<dynamicRefine2DFvMesh>(mesh_)
    );
}


Foam::geometricVofExt::SimPLIC::plicSurface
Foam::geometricVofExt::SimPLIC::reconstruction::interface()
{
    DynamicList<point> interfacePts(0.5 * mesh_.nPoints());
    DynamicList<face> interfaceFaces(0.5 * mesh_.nFaces());
    bitSet interfaceCellAddressing(mesh_.nCells());

    interfacePts.clear();
    interfaceFaces.clear();

    forAll(mixedCells_, i)
    {
        const label cellI(mixedCells_[i]);

        const vector& normalI(interfaceN_.primitiveField()[cellI]);
        const scalar distanceI(interfaceD_.primitiveField()[cellI]);
        const label cellStatusI
        (
            cutCell_.calcSubCell(cellI, normalI, distanceI, false)
        );

        if (cellStatusI == 0)
        {
            interfaceFaces.append
            (
                face
                (
                    identity
                    (
                        cutCell_.interfacePoints().size(),
                        interfacePts.size()
                    )
                )
            );

            interfacePts.append(cutCell_.interfacePoints());

            interfaceCellAddressing.set(cellI);
        }
    }

    labelList meshCells(interfaceCellAddressing.sortedToc());

    // Transfer to mesh storage
    pointField pts(std::move(interfacePts));
    faceList faces(std::move(interfaceFaces));

    return plicSurface(std::move(pts), std::move(faces), std::move(meshCells));
}


Foam::geometricVofExt::SimPLIC::reconstructedSubcellFaces
Foam::geometricVofExt::SimPLIC::reconstruction::subCellFaces()
{
    DynamicList<point> subCellPts(0.5 * mesh_.nPoints());
    DynamicList<face> subCellFaces(0.5 * mesh_.nFaces());
    bitSet interfaceCellAddressing(mesh_.nCells());

    subCellPts.clear();
    subCellFaces.clear();

    forAll(mixedCells_, i)
    {
        const label cellI(mixedCells_[i]);

        const vector& normalI(interfaceN_.primitiveField()[cellI]);
        const scalar distanceI(interfaceD_.primitiveField()[cellI]);
        const label cellStatusI
        (
            cutCell_.calcSubCell(cellI, normalI, distanceI, false)
        );

        if (cellStatusI == 0)
        {
            forAll(cutCell_.subCellFaces(), i)
            {
                face facei(cutCell_.subCellFaces()[i]);

                forAll(facei, pti)
                {
                    facei[pti] += subCellPts.size();
                }

                subCellFaces.append(facei);
            }

            subCellPts.append(cutCell_.subCellPoints());

            interfaceCellAddressing.set(cellI);
        }
    }

    labelList meshCells(interfaceCellAddressing.sortedToc());

    // Transfer to mesh storage
    pointField pts(std::move(subCellPts));
    faceList faces(std::move(subCellFaces));

    return reconstructedSubcellFaces
    (
        std::move(pts),
        std::move(faces),
        std::move(meshCells)
    );
}


// ************************************************************************* //