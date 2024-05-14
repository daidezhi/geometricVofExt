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
    createIcosphere

Description
    Create an icosphere with given centre, radius and subdivision level,
    and write to a file with an output name.

Usage
    createIcosphere [OPTIONS] <output>

    Arguments:
        <output>          The output icosphere surface file
    Options:
        -centre <point>   Centre of the icosphere (Default: (0 0 0))
        -radius <scalar>  Radius of the icosphere (Default: 1)
        -subdiv <label>   Subdivision level of the icosphere (Default: 3)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "OSspecific.H"
#include "MeshedSurfaces.H"
#include "mergePoints.H"
#include "fileName.H"
#include "foamVtkSurfaceWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create an icosphere with given centre, radius and subdivision level, "
        "and write to a file with an output name."
    );

    // Disable and add some options
    argList::noParallel();          // Disable parallel function
    argList::noFunctionObjects();   // Don't use function objects

    argList::addArgument
    (
        "output", "The output icosphere surface file"
    );

    argList::addOption
    (
        "centre",
        "point",
        "Centre of the icosphere (Default: (0 0 0))"
    );

    argList::addOption
    (
        "radius",
        "scalar",
        "Radius of the icosphere (Default: 1)"
    );

    argList::addOption
    (
        "subdiv",
        "label",
        "Subdivision level of the icosphere (Default: 3)"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    const fileName exportName(args.get<fileName>(1));

    const scalar startTime(runTime.elapsedCpuTime());

    const point centre
    (
        args.found("centre") ? args.get<point>("centre") : point::zero
    );

    const scalar radius
    (
        args.found("radius") ? args.get<scalar>("radius") : 1.0
    );

    const label subdivision
    (
        args.found("subdiv") ? args.get<label>("subdiv") : 3
    );

    Info<< nl << "Creating an icosphere ..." << nl
        << "    centre    = " << centre << nl
        << "    radius    = " << radius << nl
        << "    level     = " << subdivision << nl
        << endl;


    // Initial points and faces
    //   from PyMesh (https://github.com/PyMesh/PyMesh)
    const scalar paraR(0.5 * (1.0 + Foam::sqrt(5.0)));
    pointField icospherePoints
    (
        std::initializer_list<point>
        {
            point(  -1.0,  paraR,    0.0), point(   1.0,  paraR,    0.0),
            point(  -1.0, -paraR,    0.0), point(   1.0, -paraR,    0.0),
            point(   0.0,   -1.0,  paraR), point(   0.0,    1.0,  paraR),
            point(   0.0,   -1.0, -paraR), point(   0.0,    1.0, -paraR),
            point( paraR,    0.0,   -1.0), point( paraR,    0.0,    1.0),
            point(-paraR,    0.0,   -1.0), point(-paraR,    0.0,    1.0)
        }
    );

    faceList icosphereFaces
    (
        std::initializer_list<face>
        {
            face(triFace( 0, 11,  5)), face(triFace( 0,  5,  1)),
            face(triFace( 0,  1,  7)), face(triFace( 0,  7, 10)),
            face(triFace( 0, 10, 11)), face(triFace( 1,  5,  9)),
            face(triFace( 5, 11,  4)), face(triFace(11, 10,  2)),
            face(triFace(10,  7,  6)), face(triFace( 7,  1,  8)),
            face(triFace( 3,  9,  4)), face(triFace( 3,  4,  2)),
            face(triFace( 3,  2,  6)), face(triFace( 3,  6,  8)),
            face(triFace( 3,  8,  9)), face(triFace( 5,  4,  9)),
            face(triFace( 2,  4, 11)), face(triFace( 6,  2, 10)),
            face(triFace( 8,  6,  7)), face(triFace( 9,  8,  1))
        }
    );

    // Subdivision
    for (label subdivLevel = 1; subdivLevel <= subdivision; subdivLevel++)
    {
        // Copy previous faces and clear
        const faceList icosphereFacesCopy(icosphereFaces);
        icosphereFaces.clear();

        // Perform subdivision for each triangle
        forAll (icosphereFacesCopy, facei)
        {
            const face& fa(icosphereFacesCopy[facei]);

            /* Compute 3 extra points by spliting half on each edge
            //           P1
            //          / \
            //   newP1 *---* newP3
            //        / \ / \
            //      P2---*---P3
            //         newP2
            */
            const point& P1(icospherePoints[fa[0]]);
            const point& P2(icospherePoints[fa[1]]);
            const point& P3(icospherePoints[fa[2]]);

            const point newP1(0.5 * (P1 + P2));
            const point newP2(0.5 * (P3 + P2));
            const point newP3(0.5 * (P1 + P3));

            const label maxPtId(icospherePoints.size() - 1);

            icospherePoints.append(newP1);
            icospherePoints.append(newP2);
            icospherePoints.append(newP3);

            icosphereFaces.append(face(triFace(    fa[0], maxPtId+1, maxPtId+3)));
            icosphereFaces.append(face(triFace(maxPtId+1,     fa[1], maxPtId+2)));
            icosphereFaces.append(face(triFace(maxPtId+1, maxPtId+2, maxPtId+3)));
            icosphereFaces.append(face(triFace(maxPtId+3, maxPtId+2,     fa[2])));
        }

        // Merge duplicated points and update point labels of faces
        labelList pointToUnique;
        inplaceMergePoints<pointField>(icospherePoints, 1e-8, false, pointToUnique);

        forAll (icosphereFaces, facei)
        {
            face& fa(icosphereFaces[facei]);
            fa[0] = pointToUnique[fa[0]];
            fa[1] = pointToUnique[fa[1]];
            fa[2] = pointToUnique[fa[2]];
        }
    }

    // Scaling
    for (point& pointi: icospherePoints)
    {
        pointi = pointi * radius / mag(pointi) + centre;
    }

    Info<< "    nVertices = " << (10*pow(4,subdivision)+2) << nl
        << "    nFaces    = " << (20*pow(4,subdivision)) << nl
        << "    nEdges    = " << (30*pow(4,subdivision)) << nl
        << endl;


    // Write to file
    Info<< "Writing to "
        << word(exportName)
        //<< (exportName.has_ext("vtk") ? " " : (".vtk "))
        << "..." << endl;

//    vtk::surfaceWriter writer
//    (
//        icospherePoints,
//        icosphereFaces,
//        //vtk::formatType::LEGACY_BINARY,
//        vtk::formatType::LEGACY_ASCII,
//        runTime.globalPath()/word(exportName)
//    );

    //writer.writeGeometry();

    const meshedSurface mSurf(icospherePoints, icosphereFaces);
    mSurf.write(exportName);

    Info<< "\nExecution time: "
        << runTime.elapsedCpuTime() - startTime
        << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //