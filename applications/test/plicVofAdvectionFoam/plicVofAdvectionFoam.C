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

Application
    plicVofAdvectionFoam

Description
    Test the geometric VOF library, SimPLIC, with a solenoidal velocity
    field given by:
        Ux = cos(twoPI t / T) (2 sin^2(PI x) sin(twoPI y) sin(twoPI z)),
        Uy = cos(twoPI t / T) (-sin(twoPI x)  sin^2(PI y) sin(twoPI z)),
        Uz = cos(twoPI t / T) (-sin(twoPI x) sin(twoPI y)  sin^2(PI z)),
    where T = 6 s.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "subCycle.H"
#include "pimpleControl.H"
#include "fvcSmooth.H"

#include "solveVofEqu.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test the geometric VOF library, SimPLIC, with a solenoidal velocity"
        "field given by:\n"
        "  Ux = cos(twoPI t / T) (2 sin^2(PI x) sin(twoPI y) sin(twoPI z)),\n"
        "  Uy = cos(twoPI t / T) (-sin(twoPI x)  sin^2(PI y) sin(twoPI z)),\n"
        "  Uz = cos(twoPI t / T) (-sin(twoPI x) sin(twoPI y)  sin^2(PI z)),\n"
        "where T = 6 s."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "alphaControls.H"

    #include "plicVof.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //