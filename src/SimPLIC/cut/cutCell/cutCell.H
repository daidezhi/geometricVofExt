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
    Foam::geometricVofExt::SimPLIC::cutCell

Description
    Class designed for clipping a mesh cell with an oriented plane.

SourceFiles
    cutCell.C

\*---------------------------------------------------------------------------*/

#ifndef cutCell_H
#define cutCell_H

#include "cutFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{

/*---------------------------------------------------------------------------*\
                        Class cutCell Declaration
\*---------------------------------------------------------------------------*/

class cutCell
{
private:

    // Private data

        //- Ref to the dynamicFvMesh whose cells to be cut with a plane
        const fvMesh& mesh_;

        //- Ref to face flatness field
        const scalarField& faceFlatness_;

        //- Reference to alpha1 (VOF) field
        const volScalarField& alpha1_;

        //- Interface plane normals
        const volVectorField& interfaceN_;

        //- Interface plane signed distances
        volScalarField& interfaceD_;

        //- Interface plane centres
        volVectorField& interfaceC_;

        //- Interface plane area vectors
        volVectorField& interfaceS_;

        //- A cutFace object to reach its face cutting functionality
        cutFace cutFace_;

        //- A cell status label taking one of the values:
        //   -1: cell is fully below the interface
        //    0: cell is cut
        //   +1: cell is fully above the interface
        label cellStatus_;


        // SubCell/Faces data

            //- List of face centres for cutFaces
            DynamicList<point> cutFaceCentres_;

            //- List of face area vectors for cutFaces
            DynamicList<vector> cutFaceAreas_;

            //- Cell centre of the subcell of cellI which is "fully submerged"
            point subCellCentre_;

            //- Volume of fully submerged subcell
            scalar subCellVolume_;

            //- Volume of Fluid for cellI (subCellVolume_ / mesh_.V()[cellI])
            scalar VOF_;

            //- Points of subcell
            DynamicList<point> subCellPoints_;

            //- Faces of subcell
            DynamicList<face> subCellFaces_;

            //- Subcell points and faces is updated or not
            bool isSubCellPointsAndFacesUpdated_;


        // Interface data

            //- Cell to cut
            label cellI_;

            //- Storage for subFace edges belonging to interface
            DynamicList<DynamicList<point>> interfaceEdges_;

            //- Points constituting the cell-interface intersection
            DynamicList<point> interfacePoints_;

            //- Face centre of the interface
            point interfaceCentre_;

            //- Interface area vector (from liquid/primary to gas/secondary)
            vector interfaceArea_;

        
        // Data if splitWarpedFace == true

            // Local cell points
            pointField localPoints_;

            // Local cell faces
            faceList localFaces_;


        // Private Member Functions

            //- Calculates area and centre of the interface
            void calcInterfaceCentreAndArea();

            //- Calculates volume and centre of the cutted cell
            void calcSubCellCentreAndVolume();

            //- Get local point field and face list of cellI_. The warped faces
            //  will be split based on the value of splitWarpedFace.
            void getLocalPointFieldAndFaceList
            (
                bool splitWarpedFace = false
            );

            //- Update subCell points and faces
            void updateSubCellPointsandFaces();


public:

    // Constructors

        //- Construct from fvMesh and a surfaceScalarField
        explicit cutCell
        (
            const fvMesh& mesh,
            const scalarField& faceFlatness,
            const volScalarField& alpha1,
            const volVectorField& interfaceNormal,
            volScalarField& interfaceDistance,
            volVectorField& interfaceCentre,
            volVectorField& interfaceArea
        );


    // Member Functions

        //- Sets internal values and returns cell status
        label calcSubCell
        (
            const label cellI,
            const vector& normal,
            const scalar distance,
            bool splitWarpedFace
        );

        //- Returns subCellCentre
        const point& subCellCentre() const noexcept
        {
            return subCellCentre_;
        }

        //- Returns subCellVolume
        scalar subCellVolume() const noexcept
        {
            return subCellVolume_;
        }

        //- Returns subcell points
        const DynamicList<point>& subCellPoints()
        {
            if (!isSubCellPointsAndFacesUpdated_)
            {
                updateSubCellPointsandFaces();
            }

            return subCellPoints_;
        }

        //- Returns subcell faces
        const DynamicList<face>& subCellFaces()
        {
            if (!isSubCellPointsAndFacesUpdated_)
            {
                updateSubCellPointsandFaces();
            }

            return subCellFaces_;
        }

        //- Returns the points of the cutting PLICface
        const DynamicList<point>& interfacePoints();

        //- Returns the centre of the cutting PLICface
        const point& interfaceCentre() const noexcept
        {
            return interfaceCentre_;
        }

        //- Returns the area normal vector of the cutting PLICface
        const vector& interfaceArea() const noexcept
        {
            return interfaceArea_;
        }

        //- Returns cellStatus
        label cellStatus() const noexcept
        {
            return cellStatus_;
        }

        //- Returns volume of fluid value
        scalar volumeOfSubCell() const noexcept
        {
            return subCellVolume_;
        }

        //- Returns volume of fluid value
        scalar volumeOfFluid() const noexcept
        {
            return VOF_;
        }

        //- Find signed distance of interface plane with given cell index,
        //  fraction value and normalized orientation vector. A switch named
        //  splitWarpedFace is provided to handle cell warped faces, which are
        //  common in polyhedral meshes.
        label findSignedDistance
        (
            const label cellI,
            bool splitWarpedFace
        );

        //- Reset internal values
        void clearStorage()
        {
            cellI_ = -1;
            cellStatus_ = -1;
            cutFaceCentres_.clear();
            cutFaceAreas_.clear();
            subCellPoints_.clear();
            subCellFaces_.clear();
            isSubCellPointsAndFacesUpdated_ = false;
            interfaceEdges_.clear();
            interfacePoints_.clear();
            interfaceCentre_ = point::zero;
            interfaceArea_ = vector::zero;
            subCellCentre_ = point::zero;
            subCellVolume_ = -10;
            VOF_ = -10;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SimPLIC
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
