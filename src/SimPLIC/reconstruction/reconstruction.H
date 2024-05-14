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
    Foam::geometricVofExt::SimPLIC::reconstruction

Description
    Reconstruct the approximated interface planes.

SourceFiles
    reconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef reconstruction_H
#define reconstruction_H

#include "fvMesh.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "dynamicOversetFvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "className.H"
#include "cutFace.H"
#include "cutCell.H"
#include "plicSurface.H"
#include "reconstructedSubcellFaces.H"
#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "DynamicField.H"
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{

/*---------------------------------------------------------------------------*\
                        Class reconstruction Declaration
\*---------------------------------------------------------------------------*/

class reconstruction
:
    public IOdictionary
{
private:

    // Private data

        // External data references

            //- Ref to the mesh
            const fvMesh& mesh_;

            //- Reference to alpha1 (VOF) field
            volScalarField& alpha1_;

            //- Often used reference to alpha1 internal field
            scalarField& alpha1In_;


        // Orientation vector calculations

            //- Pointer to reconstruction scheme
            //autoPtr<orientationSchemes> orientation_;

            //- Valid orientation method
            enum validOrientationMethod_
            {
                ALPHAGRAD,
                ISOALPHAGRAD,
                ISORDF
            };

            //- Valid orientation method map
            static const std::map<Foam::word, validOrientationMethod_>
                validOrientationMethods_;

            //- Orientation vector calculation method
            word orientationMethod_;

            //- Map iterator of orientation vector calculation method
            std::map<word, validOrientationMethod_>::const_iterator
                orientationMethodMapIter_;

            //- isoRDF controls and tolerances
            scalar isoRDFTol_;
            scalar isoRDFRelTol_;
            label isoRDFIterations_;


        // Switches/tolerances

            //- Tolerance for marking mixed cells:
            //  Those with mixedCellTol_ < alpha1 < 1 - mixedCellTol_
            scalar mixedCellTol_;

            //- Write reconstruction fields or not
            bool writePlicFields_;

            //- Face flatness field
            scalarField faceFlatness_;

            //- A switch to turn on/off spliting warped faces
            bool splitWarpedFace_;

            //- Number of orientation iterations
            label nOrientationIters_;

            //- Switch to map alpha field geometrically
            bool mapAlphaField_;


        // Interface fields

            //- Interface plane normals
            volVectorField interfaceN_;

            //- Interface plane signed distances
            volScalarField interfaceD_;

            //- Interface plane centres
            volVectorField interfaceC_;

            //- Interface plane area vectors
            volVectorField interfaceS_;

            //- Face flatness
            surfaceScalarField zeta_;


        // Execution performance

            //- Time spent reconstructing the interface planes
            scalar reconstructionTime_;

            //- Time spent mapping the fraction field
            scalar alphaMappingTime_;


        // Cell and face cutting

            //- List of mixed cell labels
            DynamicList<label> mixedCells_;

            //- List of mixed cell status
            //   -1: cell is fully below the interface
            //    0: cell is cut
            //   +1: cell is fully above the interface
            DynamicList<label> cellStatus_;

            //- A cutFace object to reach its face cutting functionality
            cutFace cutFace_;

            //- A cutCell object to reach its cell cutting functionality
            cutCell cutCell_;


    // Private Member Functions

        //- No copy construct
        reconstruction(const reconstruction&) = delete;

        //- No copy assignment
        void operator=(const reconstruction&) = delete;


        // Reconstruction functions

            //- Calculate orientation vectors with alpha gradient from regular
            //  gradient schemes, which are controlled by grad(alpha1)
            //  in fvSchemes->gradSchemes
            void calcInterfaceNFromRegAlphaGrad();

            //- Calculate orientation vectors with isoAlphaGrad method (LS)
            void calcInterfaceNFromIsoAlphaGrad();


            // IsoRDF method

                //- Compute gradient at the surfaces
                void isoInterfaceGrad
                (
                    const volScalarField& phi,
                    boolList& isMixedCell,
                    DynamicField<vector>& interfaceNormal
                );

                //- Calculate orientation vectors with isoRDF method
                void calcInterfaceNFromIsoRDF();


            //- Update face flatness field
            void updateFaceFlatness();

            //- Clear out interface data
            void clearInterfaceData()
            {
                mixedCells_.clear();
                cellStatus_.clear();
            }


public:

    //- Runtime type information
    TypeName("reconstruction");


    //- Construct with given alpha field
    reconstruction
    (
        volScalarField& alpha1,
        const dictionary& dict
    );

    //- Destructor
    virtual ~reconstruction() = default;


    // Member functions

        //- Clear storage and get mixed cells
        void initialize();

        //- Reconstruct interface in mixed cells
        void reconstruct();

        //- Geometric interpolation of the alpha field in dynamic meshes
        void mapAlphaField();

        // Cell type determination functions

            //- Determine if a cell is a mixed cell
            bool isAMixedCell(const label cellI) const
            {
                return
                (
                    mixedCellTol_ < alpha1In_[cellI]
                 && alpha1In_[cellI] < 1.0 - mixedCellTol_
                );
            }

            //- Determine if a cell is an empty cell
            bool isAnEmptyCell(const label cellI) const
            {
                return
                (
                    alpha1In_[cellI] <= mixedCellTol_
                );
            }

            //- Determine if a cell is a full cell
            bool isAFullCell(const label cellI) const
            {
                return
                (
                    alpha1In_[cellI] >= 1.0 - mixedCellTol_
                );
            }


        // Access functions

            //- Return alpha field
            const volScalarField& alpha() const noexcept
            {
                return alpha1_;
            }

            //- Return mapAlphaField_ or not
            const bool willMapAlphaField() const
            {
                return mapAlphaField_;
            }

            //- Return interfaceN field
            const volVectorField& interfaceN() const
            {
                return interfaceN_;
            }

            //- Return interfaceN field
            volVectorField& interfaceN()
            {
                return interfaceN_;
            }

            //- Return interfaceD field
            const volScalarField& interfaceD() const
            {
                return interfaceD_;
            }

            //- Return interfaceD field
            volScalarField& interfaceD()
            {
                return interfaceD_;
            }

            //- Return interfaceC field
            const volVectorField& interfaceC() const
            {
                return interfaceC_;
            }

            //- Return interfaceC field
            volVectorField& interfaceC()
            {
                return interfaceC_;
            }

            //- Return interfaceS field
            const volVectorField& interfaceS() const
            {
                return interfaceS_;
            }

            //- Return interfaceS field
            volVectorField& interfaceS()
            {
                return interfaceS_;
            }

            //- Return alpha field
            const scalarField& faceFlatness() const
            {
                return faceFlatness_;
            }

            //- Return mixedCells
            const DynamicList<label>& mixedCells() const
            {
                return mixedCells_;
            }

            //- Return cellStatus
            const DynamicList<label>& cellStatus() const
            {
                return cellStatus_;
            }

            //- Return time in reconstruction step
            scalar reconstructionTime() const noexcept
            {
                return reconstructionTime_;
            }

            //- Return time in alpha mapping step
            scalar alphaMappingTime() const noexcept
            {
                return alphaMappingTime_;
            }

            //- Generated PLIC interface surface points and faces.
            //  It should be noted that the points are disconnected!
            plicSurface interface();

            //-
            reconstructedSubcellFaces subCellFaces();

            //-
            cutCell& cellCutter()
            {
                return cutCell_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SimPLIC
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //