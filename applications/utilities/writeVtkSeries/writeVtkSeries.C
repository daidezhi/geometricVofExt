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
    writeVtkSeries

Description
    Generate and write a VTK series for a `surfaces` function object
    named <funcName>. The output file is `<funcName>.vtp.series` or
    `<funcName>.vtk.series`, which depends on the VTK file extension.

    The `<funcName>.vtp.series` file format is a simple JSON format with
    the following type of content:
    \verbatim
    {
      "file-series-version" : "1.0",
      "files": [
        { "name" : "postProcessing/<funcName>/1/<faceName>.vtp", "time" : 1 },
        { "name" : "postProcessing/<funcName>/2/<faceName>.vtp", "time" : 2 },
        { "name" : "postProcessing/<funcName>/3/<faceName>.vtp", "time" : 3 },
      ]
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fileOperation.H"
#include "functionObject.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate and write a VTK series for a `surfaces` function object "
        "named <funcName>. The output file is `<funcName>.vtp.series` or "
        "`<funcName>.vtk.series`, which depends on the VTK file extension."
    );

    argList::noParallel();          // Disable parallel function
    argList::noFunctionObjects();   // Don't use function objects

    argList::addArgument
    (
        "funcName", "Name of the `surfaces` function object"
    );

    argList::addArgument
    (
        "faceName", "Name within the `surfaces` sub dictionary"
    );

    argList args(argc, argv);

    const word funcName(args.get<fileName>(1));
    const word faceName(args.get<fileName>(2));

    const fileName postDir
    (
        functionObject::outputPrefix/funcName
    );

    if (fileHandler().filePath(postDir).empty())
    {
        Info << endl;

        WarningInFunction
            << "Path " << postDir << " is not existing!"
            << nl << endl;
    }
    else
    {
        Info << nl << "Found path " << postDir << nl << endl;

        const instantList timeDirs(fileHandler().findTimes(postDir, word("")));

        if (timeDirs.size() == 0)
        {
            WarningInFunction
                << "Path " << postDir << " is empty!"
                << nl << endl;
        }
        else
        {
            DynamicList<word> timeNames;
            DynamicList<word> fileNames;

            for (auto& timeDirI: timeDirs)
            {
                for (fileName& fileI: fileHandler().readDir(postDir/timeDirI.name()))
                {
                    if (fileI.lessExt() == faceName and (fileI.ext() == "vtp" or fileI.ext() == "vtk"))
                    {
                        timeNames.append(timeDirI.name());

                        fileNames.append(postDir/timeDirI.name()/fileI);

                        break;
                    }
                }
            }

            Info << "Found " << fileNames.size() << " instances" << nl << endl;

            const word seriesName
            (
                funcName + "." + faceName + "." + fileName(fileNames[0].ext()) +".series"
            );

            Info << "Write VTK series to '" << seriesName << "'" << nl << endl;

            autoPtr<OFstream> osPtr(autoPtr<OFstream>::New(seriesName));
            Ostream& os(osPtr());

            // Begin file-series (JSON)
            os << "{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n";

            forAll(timeNames, timei)
            {
                os  << "    { \"name\" : \"" << fileNames[timei]
                    << "\", \"time\" : " << timeNames[timei] << " }";

                os  << ',' << nl;
            }

            os << "  ]\n}\n";
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
