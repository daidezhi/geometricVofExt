/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "reconstructedSubcellFaces.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricVofExt::SimPLIC::reconstructedSubcellFaces::reconstructedSubcellFaces
(
    const pointField& pointLst,
    const faceList& faceLst,
    const labelList& meshCells
)
:
    meshedSurface(pointLst, faceLst, surfZoneList()),
    meshCells_(meshCells)
{}


Foam::geometricVofExt::SimPLIC::reconstructedSubcellFaces::reconstructedSubcellFaces
(
    pointField&& pointLst,
    faceList&& faceLst,
    labelList&& meshCells
)
:
    meshedSurface
    (
        std::move(pointLst),
        std::move(faceLst),
        surfZoneList()
    ),
    meshCells_(std::move(meshCells))
{}


// ************************************************************************* //