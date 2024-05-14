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

Typedefs
    CGALKernel

Description
    CGAL data structures used by geometricVofExt.

    Use Exact_predicates_inexact_constructions kernel.

\*---------------------------------------------------------------------------*/

#ifndef CGALKernel_H
#define CGALKernel_H

// Include uint.H before CGAL headers to define __STDC_LIMIT_MACROS
#include "uint.H"

/*
#ifdef CGAL_INEXACT
    // Fast kernel using a double as the storage type
    #include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
    typedef CGAL::Exact_predicates_inexact_constructions_kernel CK;
#else
    // Very robust but expensive kernel
    #include "CGAL/Exact_predicates_exact_constructions_kernel.h"
    typedef CGAL::Exact_predicates_exact_constructions_kernel CK;
#endif
*/

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
//#include "CGAL/Exact_predicates_exact_constructions_kernel.h"
#include "CGAL/Plane_3.h"
#include "CGAL/Surface_mesh.h"
//#include "CGAL/boost/graph/graph_traits_Surface_mesh.h"
#include "CGAL/Polygon_mesh_processing/corefinement.h"
//#include "CGAL/Polygon_mesh_processing/stitch_borders.h"
#include "CGAL/Polygon_mesh_processing/measure.h"
#include "CGAL/Polyhedron_3.h"
#include "CGAL/Polyhedron_incremental_builder_3.h"
#include "CGAL/Side_of_triangle_mesh.h"
#include "CGAL/Polygon_mesh_processing/clip.h"
#include "CGAL/Polygon_mesh_processing/repair.h"
#include "CGAL/Polygon_mesh_processing/intersection.h"
#include "CGAL/Bbox_3.h"
#include "CGAL/Polygon_mesh_processing/bbox.h"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef CGAL::Exact_predicates_inexact_constructions_kernel CK;
typedef CK::Point_3                                         CPoint;
typedef CK::Plane_3                                         CPlane;
typedef CGAL::Surface_mesh<CPoint>                          CMesh;
typedef CMesh::Vertex_index                                 CMVIndex;
typedef CMesh::Face_iterator                                CMFIter;
typedef CMesh::Vertex_around_face_circulator                CMVFCirt;
typedef CMesh::size_type                                    CMSzType;
typedef CGAL::Polyhedron_3<CK>                              CPolyhedron;
typedef CGAL::Side_of_triangle_mesh<CMesh, CK>              CSideOfMesh;

namespace PMP = CGAL::Polygon_mesh_processing;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //