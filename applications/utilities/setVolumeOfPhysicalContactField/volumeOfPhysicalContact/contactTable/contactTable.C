/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Dezhi Dai, Argonne National Laboratory (ANL)
-------------------------------------------------------------------------------
License
    This file is part of to OpenFOAM.

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


#include "contactTable.H"

#include <vector>
#include <string>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::geometricVofExt::contactTable::typeName = "contactTable";


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::geometricVofExt::contactTable::contactTable()
:
    contactFile_(""),
    nContacts_(0),
    contacts_(dictionary())
{}


Foam::geometricVofExt::contactTable::contactTable
(
    const Foam::IOdictionary& setVolumeOfPhysicalContactField
)
:
    contactFile_(setVolumeOfPhysicalContactField.get<fileName>("contactFile")),
    nContacts_(0),
    contacts_(dictionary())
{
    csv2::Reader<> contactFileReader;

    if (contactFileReader.mmap(contactFile_.expand().c_str()))
    {
        nContacts_ = label(contactFileReader.rows()) - 1;

        auto header = contactFileReader.header();

        label i = 0;
        for (const auto& rowI: contactFileReader)
        {
            if (i >= nContacts_) break;

            dictionary dictContactI;

            std::vector<std::string> cellsInRowI;
            for (const auto& cellI: rowI)
            {
                std::string val;
                cellI.read_value(val);
                cellsInRowI.push_back(val);
            }

            dictContactI.add("pbLeft",     std::stoi(cellsInRowI[0]));
            dictContactI.add("pbRight",    std::stoi(cellsInRowI[1]));
            dictContactI.add("pbDiameter", std::stod(cellsInRowI[2]));
            dictContactI.add
            (
                "pbRadius", 0.5 * dictContactI.get<scalar>("pbDiameter")
            );

            dictContactI.add
            (
                "point1",
                Foam::point
                (
                    std::stod(cellsInRowI[3]),
                    std::stod(cellsInRowI[4]),
                    std::stod(cellsInRowI[5]))
            );

            dictContactI.add
            (
                "point2",
                Foam::point
                (
                    std::stod(cellsInRowI[6]),
                    std::stod(cellsInRowI[7]),
                    std::stod(cellsInRowI[8])
                )
            );

            dictContactI.add("chamferSize", std::stod(cellsInRowI[9]));
            dictContactI.add("contactSize", std::stod(cellsInRowI[10]));

            dictContactI.add
            (
                "radius", // physical contact radius, used for creating cylinder
                dictContactI.get<scalar>("pbRadius") * dictContactI.get<scalar>("contactSize")
            );

            dictContactI.add("type",     "cylinder");
            dictContactI.add("vertices", label(32));
            dictContactI.add("mode",     "add");

            contacts_.add(word("contact_"+std::to_string(i)), dictContactI);

            i++;
        }
    }
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //


// ************************************************************************* //