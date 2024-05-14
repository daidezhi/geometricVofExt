/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Dezhi Dai
-------------------------------------------------------------------------------
License
    This file is part of plicVof which is an extension to OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "isoSurfaceVolume.H"
#include "volPointInterpolation.H"
#include "cutCellIso.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::isoSurfaceVolume::isoSurfaceVolume
(
    const volScalarField& alpha1,
    const List<scalar>& isoValues
)
:
    alpha1_(alpha1),
    isoValues_(isoValues),
    volumes_(isoValues.size(), 0)
{
    volPointInterpolation vpi(alpha1_.mesh());

    pointScalarField alphap(vpi.interpolate(alpha1_));

    cutCellIso isoCellCutter(alpha1_.mesh(), alphap.ref());

    const scalarField& alphaIn(alpha1_.primitiveField());

    forAll(isoValues_, valuei)
    {
        volumes_[valuei] = 0.0;

        forAll(alphaIn, celli)
        {
            const label celliStatus
            (
                isoCellCutter.calcSubCell(celli, isoValues_[valuei])
            );

            if (celliStatus != 1)
            {
                volumes_[valuei] += isoCellCutter.subCellVolume();
            }
        }
    }
}


// ************************************************************************* //
