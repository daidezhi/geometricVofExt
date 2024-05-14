/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::geometricVofExt::SimPLIC::plicSurface

Description
    Mesh surface of PLIC faces for post-processing.

SourceFiles
    plicSurface.C

\*---------------------------------------------------------------------------*/

#ifndef plicSurface_H
#define plicSurface_H

#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geometricVofExt
{
namespace SimPLIC
{

/*---------------------------------------------------------------------------*\
                        Class plicSurface Declaration
\*---------------------------------------------------------------------------*/

class plicSurface
:
    public meshedSurface
{
private:

    // Private data

        //- For every face the original cell in mesh
        labelList meshCells_;


public:

    // Constructors

        //- 
        plicSurface
        (
            const pointField& pointLst,
            const faceList& faceLst,
            const labelList& meshCells
        );

        //-
        plicSurface
        (
            pointField&& pointLst,
            faceList&& faceLst,
            labelList&& meshCells
        );


    // Member Functions

        // Access functions

            //- For every face, the original cell in mesh
            const labelList& meshCells() const
            {
                return meshCells_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SimPLIC
} // End namespace geometricVofExt
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //