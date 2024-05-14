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
    Foam::geometricVofExt::CGALVof::CGALSurfaceMesh

Description
    A CGAL surface mesh.

SourceFiles
    CGALSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef CGALSurfaceMesh_H
#define CGALSurfaceMesh_H

#include "CGALKernel.H"
#include "dynamicFvMesh.H"
//#include "volFields.H"
//#include "shapeMesh.H"
#include "MeshedSurfaces.H"
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace CGALVof
{

/*---------------------------------------------------------------------------*\
                        class CGALSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class CGALSurfaceMesh
{
private:

    // Private data

        //- CGAL polyhedron
        CMesh surface_;

        //- Shape enumerations
        enum validShape_
        {
            PLANE,
            BOX,
            ROTATEDBOX,
            SPHERE,
            CYLINDER,
            CONE,
            TORUS,
            FILE
        };

        //- Valid Shape map
        static const std::map<Foam::word, validShape_> validShapes_;


    // Private Member Functions

        //- Create a box
        void createBox(const dictionary& shape);

        //- Create an icosphere
        void createSphere(const dictionary& shape);

        //- Create a cylinder
        void createCylinder(const dictionary& shape);

        //- Create a cone
        void createCone(const dictionary& shape);

        //- Create from file
        void createFromFile(const dictionary& shape);


public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Null contructor
        CGALSurfaceMesh();

        //- Construct from file name (uses extension to determine type)
        //CGALSurfaceMesh(const fileName& fName);

        //- Construct from shape dictionary
        CGALSurfaceMesh(const dictionary& shape);

        //- Construct from a cell in a finite volume mesh
        CGALSurfaceMesh
        (
            const fvMesh& mesh,
            const label cellI,
            const word triangulationScheme = "fan"
        );

        //- Construct from given point and face list
        CGALSurfaceMesh
        (
            const DynamicList<point>& pList,
            const DynamicList<face>& fList
        );

        //- Construct from a meshedSurface
        CGALSurfaceMesh(const meshedSurface& mSurf);

        //- Construct from a CGAL surface mesh
        CGALSurfaceMesh(const CMesh& cgalSurfMesh);


    //- Destructor
    virtual ~CGALSurfaceMesh() = default;


    // Member functions

        //- Return CGAL polyhedron
        const CMesh& surface() const
        {
            return surface_;
        }

        //- Return CGAL polyhedron
        CMesh& surface()
        {
            return surface_;
        }

        //- Return volume
        scalar volume() const;

        //- Convert from OpenFOAM point field and face list
        void convert(const pointField& cPts, const faceList& fList);

        //- Return OpenFOAM point field
        pointField getPointField() const;

        //- Return OpenFOAM face list
        faceList getFaceList() const;

        //- Calculate clipped shape
        bool clip(const vector& normal, const scalar distance);

        //- Write to file
        void write(const fileName& fName) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGALVof
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
