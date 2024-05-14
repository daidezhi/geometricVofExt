/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2017 OpenCFD Ltd.
    Modified code Copyright (C) 2019-2020 DLR
    Modified code Copyright (C) 2021 Johan Roenby
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

#include "advection.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "localMin.H"

// ************************************************************************* //

template<typename Type>
Type Foam::geometricVofExt::SimPLIC::advection::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label faceI
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        return f.primitiveField()[faceI];
    }
    else
    {
        const polyBoundaryMesh& pbm(mesh_.boundaryMesh());

        // Boundary face. Find out which face of which patch
        const label patchi(pbm.patchID()[faceI - mesh_.nInternalFaces()]);

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << faceI
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp(pbm[patchi]);
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return pTraits<Type>::zero;
        }

        const label patchFacei(pp.whichFace(faceI));
        return f.boundaryField()[patchi][patchFacei];
    }
}


template<typename Type>
void Foam::geometricVofExt::SimPLIC::advection::setFaceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label faceI,
    const Type& value
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        f.primitiveFieldRef()[faceI] = value;
    }
    else
    {
        const polyBoundaryMesh& pbm(mesh_.boundaryMesh());

        // Boundary face. Find out which face of which patch
        const label patchi(pbm.patchID()[faceI - mesh_.nInternalFaces()]);

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << faceI
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp(pbm[patchi]);
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return;
        }

        const label patchFacei(pp.whichFace(faceI));

        f.boundaryFieldRef()[patchi][patchFacei] = value;
    }
}


template<class SpType, class SuType>
void Foam::geometricVofExt::SimPLIC::advection::limitFlux
(
    const SpType& Sp,
    const SuType& Su
)
{
    const scalar aTol(100.0*SMALL);                 // Tolerance
    scalar maxAlphaMinus1(gMax(alpha1_) - 1.0);     // max(alphaNew - 1);
    scalar minAlpha(gMin(alpha1_));                 // min(alphaNew);
    const label nOvershoots(20);            // sum(pos0(alphaNew - 1 - aTol));

    const labelList& owner(mesh_.faceOwner());
    const labelList& neighbour(mesh_.faceNeighbour());

    Info<< "SimPLIC::advection: Before conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    surfaceScalarField dVfCorrectionValues("dVfCorrectionValues", dVf_ * 0.0);

    bitSet needBounding(mesh_.nCells(), false);
    needBounding.set(reconstructor_.mixedCells());

    extendMarkedCells(needBounding);

    // Loop number of bounding steps
    for (label n = 0; n < nAlphaBounds_; n++)
    {
        if (maxAlphaMinus1 > aTol || minAlpha < -aTol)
        {
            DynamicList<label> correctedFaces(3 * nOvershoots);
            dVfCorrectionValues = dimensionedScalar("0", dimVolume, 0.0);
            boundFlux
            (
                needBounding,
                dVfCorrectionValues,
                correctedFaces,
                Sp,
                Su
            );

            correctedFaces.append
            (
                syncProcPatches(dVfCorrectionValues, phi_, true)
            );

            labelHashSet alreadyUpdated;
            for (const label facei : correctedFaces)
            {
                if (alreadyUpdated.insert(facei))
                {
                    checkIfOnProcPatch(facei);
                    const label own(owner[facei]);
                    scalar Vown(mesh_.V()[own]);

                    alpha1_[own] -=
                        faceValue(dVfCorrectionValues, facei) / Vown;

                    if (mesh_.isInternalFace(facei))
                    {
                        const label nei(neighbour[facei]);
                        scalar Vnei(mesh_.V()[nei]);

                        alpha1_[nei] +=
                            faceValue(dVfCorrectionValues, facei) / Vnei;
                    }

                    // Change to treat boundaries consistently
                    const scalar corrVf =
                        faceValue(dVf_, facei)
                      + faceValue(dVfCorrectionValues, facei);

                    setFaceValue(dVf_, facei, corrVf);
                }
            }
            syncProcPatches(dVf_, phi_);
        }
        else
        {
            break;
        }

        // if (isA<dynamicOversetFvMesh>(mesh_))
        // {
        //     mesh_.interpolate(alpha1_);

        //     applyBruteForceBounding();
        // }

        maxAlphaMinus1 = gMax(alpha1_) - 1.0;   // max(alphaNew - 1);
        minAlpha = gMin(alpha1_);               // min(alphaNew);
    }

    Info<< "SimPLIC::advection: After  conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    alpha1_.correctBoundaryConditions();
}


template<class SpType, class SuType>
void Foam::geometricVofExt::SimPLIC::advection::boundFlux
(
    const bitSet& nextToInterface,
    surfaceScalarField& dVfCorrectionValues,
    DynamicList<label>& correctedFaces,
    const SpType& Sp,
    const SuType& Su
)
{
    const scalar aTol(100.0*SMALL);                     // Tolerance
    const scalar dt(mesh_.time().deltaTValue());
    const scalar rDeltaT(1.0 / dt);

    correctedFaces.clear();

    const scalarField& meshV(mesh_.cellVolumes());
    const volScalarField& alphaOld(alpha1_.oldTime());

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());

    // Loop through alpha cell centred field
    for (label celli: nextToInterface)
    {
        if (alpha1_[celli] < -aTol || alpha1_[celli] > 1.0 + aTol)
        {
            scalar Vi(meshV[celli]);

            scalar alphaOvershoot =
                pos0(alpha1_[celli] - 1.0) * (alpha1_[celli] - 1.0)
              + neg0(alpha1_[celli]) * alpha1_[celli];

            scalar fluidToPassOn(alphaOvershoot * Vi);
            label nFacesToPassFluidThrough(1);

            bool firstLoop(true);

            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            for (label iter=0; iter<10; iter++)
            {
                if (mag(alphaOvershoot) < aTol || nFacesToPassFluidThrough == 0)
                {
                    break;
                }

                facesToPassFluidThrough.clear();
                dVfmax.clear();
                phi.clear();

                // Find potential neighbour cells to pass surplus phase to
                setDownwindFaces(celli, downwindFaces);

                scalar dVftot(0);
                nFacesToPassFluidThrough = 0;

                for (const label facei : downwindFaces)
                {
                    const scalar phif(faceValue(phi_, facei));

                    const scalar dVff =
                        faceValue(dVf_, facei)
                      + faceValue(dVfCorrectionValues, facei);

                    const scalar maxExtraFaceFluidTrans =
                        mag(pos0(fluidToPassOn) * phif * dt - dVff);

                    // dVf has same sign as phi and so if phi>0 we have
                    // mag(phi_[facei]*dt) - mag(dVf[facei]) = phi_[facei]*dt
                    // - dVf[facei]
                    // If phi < 0 we have mag(phi_[facei]*dt) -
                    // mag(dVf[facei]) = -phi_[facei]*dt - (-dVf[facei]) > 0
                    // since mag(dVf) < phi*dt
                    if (maxExtraFaceFluidTrans / Vi > aTol)
                    {
                        facesToPassFluidThrough.append(facei);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif * dt);
                    }
                }

                forAll(facesToPassFluidThrough, fi)
                {
                    const label facei(facesToPassFluidThrough[fi]);
                    scalar fluidToPassThroughFace =
                        mag(fluidToPassOn) * mag(phi[fi] * dt) / dVftot;

                    nFacesToPassFluidThrough +=
                        pos0(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(dVfCorrectionValues, facei);

                    dVff +=
                        sign(phi[fi]) * sign(fluidToPassOn)
                      * fluidToPassThroughFace;

                    setFaceValue(dVfCorrectionValues, facei, dVff);

                    if (firstLoop)
                    {
                        checkIfOnProcPatch(facei);
                        correctedFaces.append(facei);
                    }
                }

                firstLoop = false;

                scalar alpha1New
                (
                    (
                        alphaOld[celli] * rDeltaT + Su[celli]
                      - netFlux(celli, dVf_) / Vi * rDeltaT
                      - netFlux(celli, dVfCorrectionValues) / Vi * rDeltaT
                    ) / (rDeltaT - Sp[celli])
                );

                alphaOvershoot =
                    pos0(alpha1New - 1.0) * (alpha1New - 1.0)
                  + neg0(alpha1New) * alpha1New;

                fluidToPassOn = alphaOvershoot * Vi;
            }
        }
    }
}


template<class SpType, class SuType>
void Foam::geometricVofExt::SimPLIC::advection::advect
(
    const SpType& Sp,
    const SuType& Su
)
{
    if (mesh_.topoChanging())
    {
        setProcessorPatches();
    }

    scalar startTime(mesh_.time().elapsedCpuTime());

    const scalar rDeltaT(1.0 / mesh_.time().deltaTValue());

    // Initialising dVf with upwind values
    // i.e. phi[facei]*alpha1[upwindCell[facei]]*dt
    // Make sure phi_ *= faceMask if overset mesh is used
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_) * mesh_.time().deltaT();

    // Calculate volumetric face transport during dt
    timeIntegratedFlux();

    // Adjust alpha for mesh motion
    if (mesh_.moving())
    {
        alpha1In_ *= (mesh_.Vsc0() / mesh_.Vsc());
    }

    // Overset
    if (isA<dynamicOversetFvMesh>(mesh_))
    {
        const volScalarField& cellMask_
        (
            mesh_.lookupObject<volScalarField>("cellMask")
        );

        surfaceScalarField faceMask_
        (
            localMin<scalar>(mesh_).interpolate(cellMask_)
        );

        dVf_ *= faceMask_;
    }

    // Advect the free surface
    alpha1_.primitiveFieldRef() =
    (
        alpha1_.oldTime().primitiveField() * rDeltaT
      + Su.field()
      - fvc::surfaceIntegrate(dVf_)().primitiveField() * rDeltaT
    ) / (rDeltaT - Sp.field());

    alpha1_.correctBoundaryConditions();

    // Adjust dVf for unbounded cells
    limitFlux(Sp, Su);

    // Apply non-conservative bounding mechanisms (clipping and snapping)
    // Note: We should be able to write out alpha before this is done!
    applyBruteForceBounding();

    advectionTime_ += (mesh_.time().elapsedCpuTime() - startTime);

    alphaPhi_ = dVf_ / mesh_.time().deltaT();
}


// ************************************************************************* //