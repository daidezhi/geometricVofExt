/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::geometricVofExt::SimPLIC::solveVofEqu::typeName = "solveVofEqu";


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::geometricVofExt::SimPLIC::solveVofEqu::solveVofEqu
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    // External data references
    mesh_(dynamic_cast<const fvMesh&>(alpha1.mesh())),
    alpha1_(alpha1),
    alpha1In_(alpha1.ref()),
    phi_(phi),
    U_(U),
    dict_(mesh_.solverDict(alpha1.name())),

    // Reconstructor
    reconstructor_(alpha1, dict_),

    // Advector
    advector_(alpha1, phi, U, reconstructor_, dict_)
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
    }
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

void Foam::geometricVofExt::SimPLIC::solveVofEqu::reconstruct()
{
    reconstructor_.reconstruct();
}


void Foam::geometricVofExt::SimPLIC::solveVofEqu::mapAlphaField()
{
    reconstructor_.mapAlphaField();
}


// ************************************************************************* //