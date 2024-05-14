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

Class
    Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation

Description
    A boolean operation processor of two surface meshes.

SourceFiles
    surfaceMeshBooleanOperation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceMeshBooleanOperation_H
#define surfaceMeshBooleanOperation_H

#include "CGALSurfaceMesh.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace CGALVof
{

/*---------------------------------------------------------------------------*\
                class surfaceMeshBooleanOperation Declaration
\*---------------------------------------------------------------------------*/

class surfaceMeshBooleanOperation
{
private:


protected:

    // Protected data

        //- Reference to the triangle surfaces
        const CGALSurfaceMesh& surface1_;
        const CGALSurfaceMesh& surface2_;

        //- Boolean action
        word action_;

        //- The result mesh of the Boolean operation
        CMesh boolSurface_;

        //- Is the Boolean operation valid
        bool isValid_;


public:

    // Static data members

        static const char* const typeName;


    // Data types

        //- Enumeration listing the possible Boolean operator types
        enum booleanOpType
        {
            UNION,          // Union of volumes
            INTERSECTION,   // Intersection of volumes
            DIFFERENCE,     // Difference of volumes
        };

        //- Valid actions
        static const Enum<booleanOpType> validActions;


    // Constructors

        //- Null contructor
        surfaceMeshBooleanOperation();

        //- Construct from two triangle surfaces
        surfaceMeshBooleanOperation
        (
            const CGALSurfaceMesh& surface1,
            const CGALSurfaceMesh& surface2,
            const word& action
        );


    // Destructor

        virtual ~surfaceMeshBooleanOperation()
        {}


    // Member functions

        //- Return is the Boolean operation valid
        const bool& isValid() const
        {
            return isValid_;
        }

        //- Compute specific Boolean operation
        void compBooleanOp();

        //- Return result mesh
        const CMesh& result() const
        {
            return boolSurface_;
        }

        //- Return volume of result surface mesh of the Boolean operation
        scalar volume() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGALVof
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //