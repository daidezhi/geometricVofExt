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
    calcSurfaceMeshNonUniformFlow

Description
    Calculate point locations of an input surface mesh at different time
    instances with the solenoidal velocity field given by:
        Ux = cos(twoPI t / T) (2 sin^2(PI x) sin(twoPI y) sin(twoPI z)),
        Uy = cos(twoPI t / T) (-sin(twoPI x)  sin^2(PI y) sin(twoPI z)),
        Uz = cos(twoPI t / T) (-sin(twoPI x) sin(twoPI y)  sin^2(PI z)),
    where T = 6 s. The default ODE solver is 4/5th Order Dormand-Prince
    Runge-Kutta method (RKDP45).

Usage
    calcSurfaceMeshNonUniformFlow [OPTIONS] <input>

    Arguments:
        <input>           The input icosphere surface file
    Options:
        -ODESolver <word> ODE solver
        -dt <scalar>      Time step

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"

#include "unitConversion.H"
#include "OSspecific.H"
#include "MeshedSurfaces.H"
#include "mergePoints.H"
#include "fileName.H"
#include "motionODE.H"
#include "foamVtkSurfaceWriter.H"
#include "writeFile.H"

using namespace Foam;

scalar calcMeshedSurfaceVol(const meshedSurface& mSurf)
{
    return gSum(mSurf.Cf() & mSurf.Sf()) / 3.0;
}


scalar calcMeshedSurfaceArea(const meshedSurface& mSurf)
{
    return gSum(mSurf.magSf());
}


scalarField calcAspectRatios(const meshedSurface& mSurf)
{
    const scalarField& magSf(mSurf.magSf());
    const pointField& points(mSurf.points());
    const faceList& faces(mSurf.surfFaces());

    scalarField aspectRatios(faces.size());

    forAll (faces, facei)
    {
        const face& fa(faces[facei]);

        scalarList edgeLengths(3);
        edgeLengths[0] = mag(points[fa[1]] - points[fa[0]]);
        edgeLengths[1] = mag(points[fa[2]] - points[fa[1]]);
        edgeLengths[2] = mag(points[fa[0]] - points[fa[2]]);

        aspectRatios[facei]
            = max(edgeLengths)*sum(edgeLengths) / (Foam::sqrt(48.0)*magSf[facei]);
    }

    return aspectRatios;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculate point locations of an input surface mesh at different time "
        "instances with the solenoidal velocity field given by:\n"
        "  Ux = cos(twoPI t / T) (2 sin^2(PI x) sin(twoPI y) sin(twoPI z)),\n"
        "  Uy = cos(twoPI t / T) (-sin(twoPI x)  sin^2(PI y) sin(twoPI z)),\n"
        "  Uz = cos(twoPI t / T) (-sin(twoPI x) sin(twoPI y)  sin^2(PI z)),\n"
        "where T = 6 s. The default ODE solver is 4/5th Order "
        "Dormand-Prince Runge-Kutta method (RKDP45)."
    );

    // Disable and add some options
    argList::noParallel();          // Disable parallel function
    argList::noFunctionObjects();   // Don't use function objects

    argList::addArgument
    (
        "input", "The input icosphere surface file"
    );

    argList::addOption
    (
        "dt",
        "scalar",
        "Time step"
    );

    argList::addOption
    (
        "ODESolver",
        "word",
        "ODE solver"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName(), false, false);

    runTime.setTime(scalar(0), 0);
    runTime.setEndTime(scalar(3));
    runTime.setDeltaT
    (
        args.found("dt") ? args.get<scalar>("dt") : 1.0,
        false
    );

    const fileName importName(args.get<fileName>(1));

    Info << nl << "Reading " << importName << " ..." << nl << endl;

    meshedSurface mSurf0(importName);
    const pointField& mSurfPoints0(mSurf0.points());
    pointField mSurfPoints(mSurf0.points());
    const faceList mSurfFaces(mSurf0.surfFaces());

    scalar shapeVol(calcMeshedSurfaceVol(mSurf0));
    scalar shapeArea(calcMeshedSurfaceArea(mSurf0));
    scalarField shapeAspectRatios(calcAspectRatios(mSurf0));

    #include "writeToFile.H"

    #include "createMotionSolver.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime.setTime
        (
            runTime.timeOutputValue() + runTime.deltaTValue(),
            runTime.timeIndex() + 1
        );

        Info<< "Time = " << runTime.timeName() << endl;

        Info << "Solving moving distances ..." << endl;

        forAll (mSurfPoints, pointi)
        {
            sStart[0] = mSurfPoints0[pointi].x();
            sStart[1] = mSurfPoints0[pointi].y();
            sStart[2] = mSurfPoints0[pointi].z();

            scalar t(tStart);
            scalarField s(sStart);

            scalar dtEst(1e-6);
            motionSolver->solve(t, runTime.timeOutputValue(), s, dtEst);

            mSurfPoints[pointi].x() = s[0];
            mSurfPoints[pointi].y() = s[1];
            mSurfPoints[pointi].z() = s[2];
        }

        const meshedSurface mSurfI(mSurfPoints, mSurfFaces);

        const fileName writePath
        (
            outputDir / importName.nameLessExt() / runTime.timeName()
        );

        if (!exists(writePath))
        {
            mkDir(writePath);
        }

        mSurfI.write(writePath / "shape.vtk");
        mSurfI.write(writePath / "shape.vtp");
        mSurfI.write(writePath / "shape.stlb");

        shapeVol = calcMeshedSurfaceVol(mSurfI);
        shapeArea = calcMeshedSurfaceArea(mSurfI);
        shapeAspectRatios = calcAspectRatios(mSurfI);

        os << Foam::setw(2) << "  ";
        os << Foam::setf(ios_base::left) << Foam::setw(10) << runTime.timeName();
        os  << Foam::setf(ios_base::left) << Foam::setw(20)
            << gMin(shapeAspectRatios);
        os  << Foam::setf(ios_base::left) << Foam::setw(25)
            << gMax(shapeAspectRatios);
        os  << Foam::setf(ios_base::left) << Foam::setw(25)
            << gAverage(shapeAspectRatios);
        os << Foam::setf(ios_base::left) << Foam::setw(20) << shapeArea;
        os << Foam::setf(ios_base::left) << Foam::setw(17) << shapeVol << nl;

        os.flush();

        runTime.printExecutionTime(Info);
    }

    const scalar relativeDistanceErr
    (
        gSum(mag(mSurfPoints - mSurfPoints0)) / 0.15
    );
    Info<< "Relative distance error at time of 3s = "
        << relativeDistanceErr
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //