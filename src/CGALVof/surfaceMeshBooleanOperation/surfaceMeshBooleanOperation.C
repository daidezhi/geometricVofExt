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

#include "surfaceMeshBooleanOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation::typeName
=
"surfaceMeshBooleanOperation";


const Foam::Enum
<
    Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation::booleanOpType
>
Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation::validActions
{
    { INTERSECTION, "intersection"},
    { UNION,        "union"},
    { DIFFERENCE,   "difference"}
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation::surfaceMeshBooleanOperation()
:
    surface1_(CGALSurfaceMesh()),
    surface2_(CGALSurfaceMesh()),
    action_(word("intersection")),
    boolSurface_(CMesh()),
    isValid_(false)
{}


Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation::surfaceMeshBooleanOperation
(
    const CGALSurfaceMesh& surface1,
    const CGALSurfaceMesh& surface2,
    const word& action
)
:
    surface1_(surface1),
    surface2_(surface2),
    action_(action),
    boolSurface_(CMesh()),
    isValid_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation::compBooleanOp()
{
    CMesh surf1(surface1_.surface());
    CMesh surf2(surface2_.surface());

    boolSurface_.clear();

    if (validActions[action_] == INTERSECTION)
    {
        try
        {
            isValid_
            =
            PMP::corefine_and_compute_intersection(surf1, surf2, boolSurface_);
        }catch( ... )
        {
            isValid_ = false;
        }
    }
    else if (validActions[action_] == UNION)
    {
        try
        {
            isValid_
            =
            PMP::corefine_and_compute_union(surf1, surf2, boolSurface_);
        }catch( ... )
        {
            isValid_ = false;
        }
    }
    else if (validActions[action_] == DIFFERENCE)
    {
        try
        {
            isValid_
            =
            PMP::corefine_and_compute_difference(surf1, surf2, boolSurface_);
        }catch( ... )
        {
            isValid_ = false;
        }
    }
}


Foam::scalar Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation::volume() const
{
    if (isValid_)
    {
        if (CGAL::is_closed(boolSurface_))
        {
            return scalar(PMP::volume(boolSurface_));
        }
        else
        {
            return scalar(0);
        }
    }
    else
    {
        return scalar(0);
    }
}


// ************************************************************************* //
