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

\*---------------------------------------------------------------------------*/

#include "cellToSurfaceMeshLocation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::geometricVofExt::CGALVof::cellToSurfaceMeshLocation::typeName
=
"cellToSurfaceMeshLocation";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricVofExt::CGALVof::cellToSurfaceMeshLocation::cellToSurfaceMeshLocation
(
    const fvMesh& mesh,
    const CGALSurfaceMesh& surface,
    const bool extraIntersectionCheck
)
:
    mesh_(mesh),
    surface_(surface),
    pointLocations_(mesh_.nPoints(), -2),
    faceLocations_(mesh_.nFaces(), -3),
    cellLocations_(mesh_.nCells(), -2),
    submergedCells_(0),
    intersectedCells_(0),
    aridCells_(0)
{
    const pointField& points(mesh_.points());
    const faceList& faces(mesh_.faces());
    const cellList& cellFaces(mesh_.cells());

    //CGAL::Bbox_3 bbox(CGAL::Polygon_mesh_processing::bbox(surface.surface()));

    CSideOfMesh inside(surface.surface());

    #pragma omp parallel for
    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        CGAL::Bounded_side pSide(inside(CPoint(pt.x(), pt.y(), pt.z())));

        if (pSide == CGAL::ON_BOUNDED_SIDE)
        {
            pointLocations_[pointI] = -1;
        }
        else if (pSide == CGAL::ON_BOUNDARY)
        {
            pointLocations_[pointI] = 0;
        }
        else
        {
            pointLocations_[pointI] = 1;
        }
    }

    #pragma omp parallel for
    forAll(faces, faceI)
    {
        label sumBoundedPts(0), sumUnboundedPts(0), sumBoundaryPts(0);

        const face& fa(faces[faceI]);

        forAll(fa, ptI)
        {
            if (pointLocations_[fa[ptI]] == -1)
            {
                sumBoundedPts += 1;
            }
            else if (pointLocations_[fa[ptI]] == 1)
            {
                sumUnboundedPts += 1;
            }
            else
            {
                sumBoundaryPts += 1;
            }
        }

        if (sumBoundaryPts == fa.size())
        {
            faceLocations_[faceI] = -2;
        }
        else
        {
            if (sumBoundedPts > 0 && sumUnboundedPts > 0)
            {
                faceLocations_[faceI] = 0;
            }
            else
            {
                if (sumBoundedPts == 0)
                {
                    faceLocations_[faceI] = 1;
                }
                else
                {
                    faceLocations_[faceI] = -1;
                }
            }
        }
    }

    #pragma omp parallel for
    forAll(cellFaces, cellI)
    {
        label sumBoundedFas(0), sumUnboundedFas(0), sumIntersectedFas(0);

        const cell& cl(cellFaces[cellI]);

        forAll(cl, faI)
        {
            if (faceLocations_[cl[faI]] == -1)
            {
                sumBoundedFas += 1;
            }
            else if (faceLocations_[cl[faI]] == 1)
            {
                sumUnboundedFas += 1;
            }
            else if (faceLocations_[cl[faI]] == 0)
            {
                sumIntersectedFas += 1;
            }
        }

        if (sumIntersectedFas > 0)
        {
            cellLocations_[cellI] = 0;
        }
        else if (sumUnboundedFas == 0)
        {
            cellLocations_[cellI] = -1;
        }
        else
        {
            cellLocations_[cellI] = 1;
        }

        if (extraIntersectionCheck)
        {
            // If mesh size is greater than tool surface mesh thickness
            if (sumIntersectedFas == 0)
            {
                const CGALSurfaceMesh cCell(mesh, cellI, word("centroid"));
                bool isIntersecting
                (
                    CGAL::Polygon_mesh_processing::do_intersect
                    (
                        surface.surface(),
                        cCell.surface()
                    )
                );

                if (isIntersecting)
                {
                    cellLocations_[cellI] = 0;
                }
            }
        }
    }

    forAll(cellFaces, cellI)
    {
        if (cellLocations_[cellI] == 0)
        {
            intersectedCells_.append(cellI);
        }
        else if (cellLocations_[cellI] == -1)
        {
            submergedCells_.append(cellI);
        }
        else
        {
            aridCells_.append(cellI);
        }
    }

    Info<< "        "
        << intersectedCells_.size() << "/" << cellFaces.size()
        << " intersected cells are detected"
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //