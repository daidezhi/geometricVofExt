/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
    Copyright (C) 2019-2020 DLR
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
    Foam::geometricVofExt::SimPLIC::cutFace

Description
    Class designed for slicing a mesh face or a face defined by given points
    with an oriented plane.

SourceFiles
    cutFace.C

\*---------------------------------------------------------------------------*/

#ifndef cutFace_H
#define cutFace_H

#include "fvMesh.H"
#include "dynamicFvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{

/*---------------------------------------------------------------------------*\
                        Class cutFace Declaration
\*---------------------------------------------------------------------------*/

class cutFace
{
private:

    // Private data

        //- Ref to the dynamicFvMesh whose faces to be cut with a plane
        const fvMesh& mesh_;

        //- Ref to face flatness field
        const scalarField& faceFlatness_;

        //- Storage for centre of subface
        point subFaceCentre_;

        //- Storage for area vector of subface
        vector subFaceArea_;

        //- Storage for subFacePoints
        DynamicList<point> subFacePoints_;

        //- Storage for interfacePoints (two endpoints of the intersected edge)
        DynamicList<point> interfacePoints_;

        //- Storage for point distances relative to a given plane
        DynamicList<scalar> pointDistances_;

        //- A face status label taking one of the values:
        //   -1: face is fully below the interface
        //    0: face is cut
        //   +1: face is fully above the interface
        label faceStatus_;


    // Private Member Functions

        //- Calculate centre and area vector of subface
        void calcSubFaceCentreAndArea();


public:

    // Constructors

        //- Construct from fvMesh and a surfaceScalarField
        explicit cutFace
        (
            const fvMesh& mesh,
            const scalarField& faceFlatness
        );


    // Member Functions

        //- Calculate cut points along edges of faceI of the mesh
        label calcSubFace
        (
            const label faceI,
            const vector& normal,
            const scalar distance
        );

        //- Calculate cut points along edges of a face with given point field
        label calcSubFace
        (
            const pointField& fPts, // face points
            const vector& normal,
            const scalar distance
        );

        //- Calculate time integrated flux for a face
        scalar timeIntegratedFaceFlux
        (
            const label faceI,
            const vector& normal,
            const scalar distance,
            const scalar Un0,
            const scalar dt,
            const scalar phi,
            const scalar magSf
        );

        //- Calculate time integrated flux for a face
        scalar timeIntegratedArea
        (
            const pointField& fPts, // face points
            const vector& normal,
            const scalar distance,
            const scalarField& pTimes,
            const scalar Un0,
            const scalar dt,
            const scalar magSf
        );

        //- Returns centre of cutted face
        const point& subFaceCentre() const noexcept
        {
            return subFaceCentre_;
        }

        //- Returns area vector of cutted face
        const vector& subFaceArea() const noexcept
        {
            return subFaceArea_;
        }

        //- Returns the cut edge of the cutted face
        const DynamicList<point>& subFacePoints() const noexcept
        {
            return subFacePoints_;
        }

        //- Returns point of the face in sorted of cutted face
        const DynamicList<point>& interfacePoints() const noexcept
        {
            return interfacePoints_;
        }

        //- Reset internal variables
        void clearStorage()
        {
            subFaceCentre_ = point::zero;
            subFaceArea_ = vector::zero;
            subFacePoints_.clear();
            interfacePoints_.clear();
            pointDistances_.clear();
            faceStatus_ = -1;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SimPLIC
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //