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
    Foam::geometricVofExt::SimPLIC::solveVofEqu

Description
    Solve the VOF equation.

SourceFiles
    solveVofEqu.C

\*---------------------------------------------------------------------------*/

#ifndef solveVofEqu_H
#define solveVofEqu_H

#include "fvMesh.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "className.H"
#include "cutFace.H"
#include "cutCell.H"
#include "reconstruction.H"
#include "advection.H"
#include "MeshedSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{

/*---------------------------------------------------------------------------*\
                        Class solveVofEqu Declaration
\*---------------------------------------------------------------------------*/

class solveVofEqu
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

            //- Reference to flux field
            const surfaceScalarField& phi_;

            //- Reference to velocity field
            const volVectorField& U_;

            //- Dictionary for solveVofEqu controls
            const dictionary& dict_;


        //- Reconstructor
        reconstruction reconstructor_;

        //- Advector
        advection advector_;


    // Private Member Functions

        //- No copy construct
        solveVofEqu(const solveVofEqu&) = delete;

        //- No copy assignment
        void operator=(const solveVofEqu&) = delete;


public:

    // Static data members

        static const char* const typeName;


    //- Construct with given alpha, phi and velocity fields
    solveVofEqu
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U
    );

    //- Destructor
    virtual ~solveVofEqu() = default;


    // Member functions

        //- Reconstruct interface in mixed cells
        void reconstruct();

        //- Advect the free surface. Updates alpha field, taking into account
        //  multiple calls within a single time step.
        template <class SpType, class SuType>
        void advect(const SpType& Sp, const SuType& Su);

        //- Geometric interpolation of the alpha field in dynamic meshes
        void mapAlphaField();


        // Access functions

            //- Return alpha field
            const volScalarField& alpha() const noexcept
            {
                return alpha1_;
            }

            //- Return control dictionary
            const dictionary& dict() const noexcept
            {
                return dict_;
            }

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const dimensionedScalar rho1,
                const dimensionedScalar rho2
            ) const
            {
                return advector_.getRhoPhi(rho1, rho2);
            }

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const volScalarField& rho1,
                const volScalarField& rho2
            )
            {
                return advector_.getRhoPhi(rho1, rho2);
            }

            //- Return alphaPhi
            const surfaceScalarField& alphaPhi() const noexcept
            {
                return advector_.alphaPhi();
            }

            //- Return reconstructor
            const reconstruction& reconstructor() const
            {
                return reconstructor_;
            }

            //- Return reconstructor
            reconstruction& reconstructor()
            {
                return reconstructor_;
            }

            //- Return time in reconstruction step
            scalar reconstructionTime() const noexcept
            {
                return reconstructor_.reconstructionTime();
            }

            //- Return time in advection step
            scalar advectionTime() const noexcept
            {
                return advector_.advectionTime();
            }

            //- Return time in advection step
            scalar alphaMappingTime() const noexcept
            {
                return reconstructor_.alphaMappingTime();
            }

            //- Write reconstructed interface polygons and sampled data
            //  $FOAM_CASE/postProcessing/plicInterface/<timeName>/interface.vtp
            //void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SimPLIC
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solveVofEquTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //