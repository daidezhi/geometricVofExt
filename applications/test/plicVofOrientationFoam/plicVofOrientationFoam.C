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
    plicVofOrientationFoam

Description
    Test orientation calculation method.

Usage
    plicVofOrientationFoam [OPTIONS]

    Options:
        -np <number>      Number of threads used in solving
        -postProcess      Execute functionObjects only

\*---------------------------------------------------------------------------*/

#include "CGALSurfaceMesh.H"
#include "cellToSurfaceMeshLocation.H"
#include "surfaceMeshBooleanOperation.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "reconstruction.H"

#include "pimpleControl.H"
#include "IOmanip.H"
#include <cmath>
#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test orientation calculation method."
    );

    argList::noParallel();          // Disable parallel function

    argList::addOption
    (
        "np",
        "number",
        "Number of threads used in solving"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"

    #include "setOpenMP.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"

    #include "createFields.H"

    geometricVofExt::SimPLIC::reconstruction reconstructor
    (
        alpha1,
        mesh.solverDict(alpha1.name())
    );

    reconstructor.reconstruct();

    #include "calcSymmetricDifference.H"

    #include "writeData.H"

    runTime.functionObjects().execute();

    Info << "End\n" << endl;

    runTime.printExecutionTime(Info);

    return 0;
}


// ************************************************************************* //