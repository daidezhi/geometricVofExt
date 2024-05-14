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
    setVofField

Description
    Set the VOF field (alpha.*) with arbitrary shapes, featuring optional
    fraction-based dynamic mesh refinement.

\*---------------------------------------------------------------------------*/

#include "CGALSurfaceMesh.H"
#include "surfaceMeshBooleanOperation.H"
#include "cellToSurfaceMeshLocation.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "dynamicOversetFvMesh.H"
#include "reconstruction.H"
#include "isoSurfaceTopo.H"
#include "volPointInterpolation.H"
#include "OSspecific.H"
#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum validMode{ADD, SUBTRACT};

const std::map<Foam::word, validMode> validModes =
{
    {"add",      validMode::ADD},
    {"subtract", validMode::SUBTRACT}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set the VOF field (alpha.*) with arbitrary shapes, featuring "
        "optional fraction-based dynamic mesh refinement."
    );

    // Disable and add some options
    argList::noParallel();          // Disable parallel function
    argList::noFunctionObjects();   // Don't use function objects

    argList::addOption("dict", "file", "Alternative setVofFieldDict");

    argList::addOption("np", "number", "Number of OpenMP threads");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "setOpenMP.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    scalar startTime(omp_get_wtime());


    #include "readControlDict.H"
    #include "createShapesAndFields.H"

    #include "setAlphaField.H"

    // Additional iterations for dynamicRefineFvMesh
    if (isA<dynamicRefineFvMesh>(mesh))
    {
        Info << "> > > Iterations for dynamicRefineFvMesh < < <" << nl << endl;

        dictionary refineDict
        (
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            ).optionalSubDict("dynamicRefineFvMeshCoeffs")
        );

        const label maxRefinement(refineDict.get<label>("maxRefinement"));

        Info << "maxRefinement = " << maxRefinement << nl << endl;

        runTime.setTime(0, 1);

        for (label iter = 0; iter < maxRefinement; iter++)
        {
            Info << "Iteration: " << (iter+1) << endl;

            mesh.update();

            alpha1In = initialValue;
            alpha1.correctBoundaryConditions();

            #include "setAlphaField.H"
        }

        runTime.setTime(0, 0);

        mesh.write();

        Info << "> > > > > > > > > > Done < < < < < < < < < < <" << nl << endl;
    }

    // Write alpha field and interfaces
    Info << "Writing field " << fieldName << nl << endl;
    alpha1.write();

    #include "exportPlicSurface.H"
    #include "exportIsoSurface.H"


    Info<< "Execution time: "
        << omp_get_wtime() - startTime
        << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
