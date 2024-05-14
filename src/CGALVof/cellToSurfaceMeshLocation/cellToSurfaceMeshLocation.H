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
    Foam::geometricVofExt::CGALVof::cellToSurfaceMeshLocation

Description
    Locations of mesh cells relative to a surface mesh.

SourceFiles
    cellToSurfaceMeshLocation.C

\*---------------------------------------------------------------------------*/

#ifndef cellToSurfaceMeshLocation_H
#define cellToSurfaceMeshLocation_H

#include "CGALSurfaceMesh.H"
#include <omp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace CGALVof
{

/*---------------------------------------------------------------------------*\
                class cellToSurfaceMeshLocation Declaration
\*---------------------------------------------------------------------------*/

class cellToSurfaceMeshLocation
{
private:


protected:

    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to the triangle surface mesh
        const CGALSurfaceMesh& surface_;

        //- Point locations
        //      -1: bounded (submerged)
        //       0: on boundary (on the surface mesh boundary)
        //       1: unbounded (not submerged)
        //      -2: not determined
        labelList pointLocations_;

        //- Face locations
        //      -1: bounded (submerged)
        //       0: intersected
        //       1: unbounded (not submerged)
        //      -2: on boundary (on the surface mesh boundary)
        //      -3: not determined
        labelList faceLocations_;

        //- Cell locations
        //      -1: bounded (submerged)
        //       0: intersected
        //       1: unbounded (not submerged)
        labelList cellLocations_;

        //- Submerged cell list
        DynamicList<label> submergedCells_;

        //- Intersected cell list
        DynamicList<label> intersectedCells_;

        //- Arid cell list
        DynamicList<label> aridCells_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellToSurfaceMeshLocation(const cellToSurfaceMeshLocation&) = delete;

        //- Disallow default bitwise copy assignment
        void operator=(const cellToSurfaceMeshLocation&) = delete;

public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct from two triangle surfaces
        cellToSurfaceMeshLocation
        (
            const fvMesh& mesh,
            const CGALSurfaceMesh& surface,
            const bool extraIntersectionCheck
        );


    // Destructor

        virtual ~cellToSurfaceMeshLocation()
        {}


    // Member functions

        //- Return submerged cell list
        const DynamicList<label>& submergedCells() const
        {
            return submergedCells_;
        }

        //- Return intersected cell list
        const DynamicList<label>& intersectedCells() const
        {
            return intersectedCells_;
        }

        //- Return intersected cell list
        const DynamicList<label>& aridCells() const
        {
            return aridCells_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGALVof
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //