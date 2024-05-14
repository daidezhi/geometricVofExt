/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
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

#include "sampledReconstructedSubcellFaces.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{
    defineTypeNameAndDebug(sampledReconstructedSubcellFaces, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledReconstructedSubcellFaces,
        word,
        reconstructedSubcellFaces
    );
}
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::updateGeometry() const
{
    const fvMesh& fvm(static_cast<const fvMesh&>(mesh()));

    // No update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Not really being used...

    // Get sub-mesh if any
    if (!subMeshPtr_ && (-1 != mesh().cellZones().findIndex(zoneNames_)))
    {
        const label exposedPatchi =
            mesh().boundaryMesh().findPatchID(exposedPatchName_);

        bitSet cellsToSelect(mesh().cellZones().selection(zoneNames_));

        subMeshPtr_.reset
        (
            new fvMeshSubset(fvm, cellsToSelect, exposedPatchi)
        );
    }

    // Clear any stored topo
    surfPtr_.clear();

    // Clear derived data
    clearGeom();

    surfPtr_.reset
    (
        new reconstructedSubcellFaces
        (
            fvm.lookupObjectRef<SimPLIC::reconstruction>
            (
                "reconstruction"
            ).subCellFaces()
        )
    );

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::sampledReconstructedSubcellFaces
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    zoneNames_(),
    exposedPatchName_(),
    surfPtr_(nullptr),
    prevTimeIndex_(-1),
    subMeshPtr_(nullptr)
{
    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }

    if (-1 != mesh.cellZones().findIndex(zoneNames_))
    {
        dict.readIfPresent("exposedPatchName", exposedPatchName_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::needsUpdate() const
{
    const fvMesh& fvm(static_cast<const fvMesh&>(mesh()));

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::expire()
{
    surfPtr_.clear();
    subMeshPtr_.clear();

    // Clear derived data
    clearGeom();

    // Already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // Force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::geometricVofExt::SimPLIC::sampledReconstructedSubcellFaces::print(Ostream& os, int level) const
{
    os  << "SimPLIC::sampledReconstructedSubcellFaces: " << name();
}


// ************************************************************************* //