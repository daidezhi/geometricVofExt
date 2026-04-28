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
    This file is part of OpenFOAM.

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

Class
    Foam::contactTable

Description
    Class of contactTable.

SourceFiles
    contactTable.C

\*---------------------------------------------------------------------------*/

#ifndef contactTable_H
#define contactTable_H

#include "IOdictionary.H"
#include "point.H"
#include "fileName.H"

#include "csv2.hpp"
#include "indicators.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{

/*---------------------------------------------------------------------------*\
                        Class contactTable Declaration
\*---------------------------------------------------------------------------*/

class contactTable
{
private:

    // Private data

        //- Contact file name
        fileName contactFile_;

        //- Number of contacts
        label nContacts_;

        //- Contact list
        dictionary contacts_;



    // Private Member Functions

        //-


public:

     // Static data members

        static const char* const typeName;


     // Constructors

        //- Null contructor
        contactTable();

        //- Construct from shape dictionary
        contactTable(const IOdictionary& setVolumeOfPhysicalContactField);


    //- Destructor
    virtual ~contactTable() = default;


    // Member functions

        //- Return number of contacts
        const label nContacts() const
        {
            return nContacts_;
        }

        //- Return contact list
        const dictionary& contacts() const
        {
            return contacts_;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //