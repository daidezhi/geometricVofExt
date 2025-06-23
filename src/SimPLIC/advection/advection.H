/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Modified code Copyright (C) 2016-2019 OpenCFD Ltd.
    Modified code Copyright (C) 2019 DLR
    Modified code Copyright (C) 2018, 2021 Johan Roenby
    Modified code Copyright (C) 2024 Dezhi Dai, ANL
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
    Foam::geometricVofExt::SimPLIC::advection

Description
    Calculates the new VOF (alpha) field after time step dt with given initial
    VOF field, a velocity field U and face fluxes phi. The fluid transport
    calculation is based on PLIC and advecting the interface across the mesh
    faces with the velocity field interpolated to the interface planes.

    This implementation is based on the isoAdvector geometric VOF method
    developed by Dr. Johan Roenby and collaborators. SimPLIC builds directly
    on the extensible structure of isoAdvector and introduces two key
    modifications:
    - An alternative interface location method based on analytical
      reconstruction algorithms for unstructured meshes. See:
      geometricVofExt::SimPLIC::cutCell::findSignedDistance().
    - A flux integration scheme using Simpson’s rule instead of explicit
      cell-face intersection tracking. See:
      geometricVofExt::SimPLIC::cutFace::timeIntegratedArea().

    References:

    Main isoAdvector idea:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    Calculation of rhoPhi:
    \verbatim
        Roenby, J., Bredmose, H., & Jasak, H. (2019).
        IsoAdvector: Geometric VOF on general meshes.
        OpenFOAM® (pp. 281-296). Springer, Cham.
    \endverbatim

    Extension to porous media flows:
    \verbatim
        Missios, K., Jacobsen, N. G., Moeller, K., & Roenby, J.
        Using the isoAdvector Geometric VOF Method for Interfacial Flows
        Through Porous Media. MARINE 2021.
    \endverbatim

    Original code supplied by Johan Roenby, DHI (2016)
    Modified Henning Scheufler, DLR

SourceFiles
    advection.C
    advectionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef advection_H
#define advection_H

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
#include "reconstruction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{

/*---------------------------------------------------------------------------*\
                        Class advection Declaration
\*---------------------------------------------------------------------------*/

class advection
{
private:

    // Private data

        // External data references

            //- Ref to the mesh
            const dynamicFvMesh& mesh_;

            //- Reference to alpha1 (VOF) field
            volScalarField& alpha1_;

            //- Often used reference to alpha1 internal field
            scalarField& alpha1In_;

            //- Reference to flux field
            const surfaceScalarField& phi_;

            //- Reference to velocity field
            const volVectorField& U_;

            //- Reference to reconstructor
            const reconstruction& reconstructor_;


        // Switches/tolerances

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Snap tolerance
            scalar snapAlphaTol_;

            //- Clip or not
            bool clip_;


        // Face cutting

            //- A cutFace object to reach its face cutting functionality
            cutFace cutFace_;

        
        // Surface liquid volume and flux

            //- Face volumetric liquid transport
            surfaceScalarField dVf_;

            //- Fraction flux
            surfaceScalarField alphaPhi_;


        // Execution performance

            //- Time spent performing interface advection
            scalar advectionTime_;


        // Boundary storages

            //- Storage for boundary faces downwind to a mixed cell
            DynamicList<label> bsFaces_;

            //- Storage for interface normals of boundary faces
            DynamicList<vector> bsn0_;

            //- Storage for interface signed distances of boundary faces
            DynamicList<scalar> bsD0_;

            //- Storage for interface velocity magnitudes of boundary faces
            DynamicList<scalar> bsUn0_;


        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicList<label> procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicList<label>> surfaceCellFacesOnProcPatches_;

    // Private Member Functions

        //- No copy construct
        advection(const advection&) = delete;

        //- No copy assignment
        void operator=(const advection&) = delete;


        // Advection functions

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(bitSet& markedCell) const;

            //- For each face calculate volumetric face transport during dt
            void timeIntegratedFlux();

            //- For a given cell return labels of fluxing-out faces 
            void setDownwindFaces
            (
                const label cellI,
                DynamicList<label>& downwindFaces
            ) const;

            //- Limit flux
            template <class SpType, class SuType>
            void limitFlux(const SpType& Sp, const SuType& Su);

            //- Bound flux
            template < class SpType, class SuType >
            void boundFlux
            (
                const bitSet& nextToInterface,
                surfaceScalarField& dVfCorrectionValues,
                DynamicList<label>& correctedFaces,
                const SpType& Sp,
                const SuType& Su
            );

            //- Given the face volume transport dVf calculates the total volume
            //  leaving a given cell. Note: cannot use dVf member because
            //  netFlux is called also for corrected dVf
            scalar netFlux
            (
                const label cellI,
                const surfaceScalarField& dVf
            ) const;

            //- Apply the bounding based on user inputs
            void applyBruteForceBounding();

            //- Clear out interface data
            void clearBoundaryData()
            {
                bsFaces_.clear();
                bsn0_.clear();
                bsD0_.clear();
                bsUn0_.clear();
            }


        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label faceI
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label faceI,
                const Type& value
            ) const;


        // Parallel run handling functions

            //- Synchronize dVf across processor boundaries using upwind value
            DynamicList<label> syncProcPatches
            (
                surfaceScalarField& dVf,
                const surfaceScalarField& phi,
                bool returnSyncedFaces = false
            );

            //- Check if the face is on processor patch and append it to the
            //  list of surface cell faces on processor patches
            void checkIfOnProcPatch(const label faceI);

            //- Set processor patches
            void setProcessorPatches();


public:

    // Static data members

        static const char* const typeName;


    //- Construct with given alpha, phi and velocity fields
    advection
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const reconstruction& reconstructor,
        const dictionary& dict
    );

    //- Destructor
    virtual ~advection() = default;


    // Member functions

        //- Advect the free surface. Updates alpha field, taking into account
        //  multiple calls within a single time step.
        template <class SpType, class SuType>
        void advect(const SpType& Sp, const SuType& Su);


        // Access functions

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const dimensionedScalar rho1,
                const dimensionedScalar rho2
            ) const
            {
                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        "rhoPhi",
                        (rho1 - rho2) * alphaPhi_ + rho2 * phi_
                    )
                );
            }

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const volScalarField& rho1,
                const volScalarField& rho2
            )
            {
                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        "rhoPhi",
                        fvc::interpolate(rho1 - rho2) * alphaPhi_
                      + fvc::interpolate(rho2) * phi_
                    )
                );
            }

            //- Return alphaPhi
            const surfaceScalarField& alphaPhi() const noexcept
            {
                return alphaPhi_;
            }

            //- Return time in reconstruction step
            scalar advectionTime() const noexcept
            {
                return advectionTime_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SimPLIC
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "advectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //