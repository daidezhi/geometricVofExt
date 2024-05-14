/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
      isoAdvector | Copyright (C) 2016-2017 DHI
    Modified code | Copyright (C) 2018-2019 Johan Roenby
    Modified code | Copyright (C) 2019-2020 DLR
    Modified code | Copyright (C) 2019-2022 Dezhi Dai
-------------------------------------------------------------------------------
License
    This file is part of plicVof, which is an extension to OpenFOAM.

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
    Foam::plicVof::isoSurfaceVolume

Description
    Calculates the new VOF (alpha) field after time step dt with given initial
    VOF field, a velocity field U and face fluxes phi. The fluid transport
    calculation is based on PLIC and advecting the interface across the mesh
    faces with the velocity field interpolated to the interface planes.

SourceFiles
    isoSurfaceVolume.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceVolume_H
#define isoSurfaceVolume_H

#include "fvMesh.H"
#include "dynamicFvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class isoSurfaceVolume Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceVolume
{
private:

    // Private data

        //- Reference to alpha1 (VOF) field
        const volScalarField& alpha1_;
        
        //- iso values
        List<scalar> isoValues_;

        //- isosurface volumes with various iso values
        List<scalar> volumes_;


public:

    // Constructors

        //- Construct from a field and iso value list
        isoSurfaceVolume
        (
            const volScalarField& alpha1,
            const List<scalar>& isoValues
        );


    // Member Functions

        // Access functions

            //- Return isoValues
            const List<scalar>& isoValues() const noexcept
            {
                return isoValues_;
            }

            //- Return isosurface volumes
            const List<scalar>& volumes() const noexcept
            {
                return volumes_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //