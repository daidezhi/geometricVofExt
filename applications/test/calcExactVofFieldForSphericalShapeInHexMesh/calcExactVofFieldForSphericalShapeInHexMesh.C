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
    calcExactVofFieldForSphericalShapeInHexMesh

Description
    Calculate the exact VOF field for a spherical shape in a hexahedral mesh.

    The exact calculation of the overlap volume between a sphere and a
    hexahedral cell is performed using the `overlap` library at
    https://github.com/severinstrobl/overlap.

    The `./Eigen` subdirectory is from `Eigen 3.4.0` at
    https://gitlab.com/libeigen/eigen/-/releases/3.4.0.

\*---------------------------------------------------------------------------*/

#include "overlap.hpp"

#include "fvCFD.H"
#include "dynamicFvMesh.H"

#include "functions.H"

#include <omp.h>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculate the exact VOF field for a spherical shape in a "
        "hexahedral mesh."
    );

    #include "addOptions.H"

    #include "setRootCase.H"

    #include "setOpenMP.H"

    #include "createTime.H"

    #include "createDynamicFvMesh.H"
    #include "checkMesh.H"

    scalar startTime(omp_get_wtime());

    #include "updateOptions.H"

    #include "createFields.H"

    #include "updateCellLocations.H"
    #include "setAlphaField.H"

    // Write alpha field and interfaces
    Info << "Writing field " << fieldName << nl << endl;
    alpha1.write();

    #include "computeShapeError.H"

    Info<< "Execution time: "
        << omp_get_wtime() - startTime
        << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //