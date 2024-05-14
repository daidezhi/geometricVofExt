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

#include "CGALSurfaceMesh.H"
#include "MeshedSurfaces.H"
#include "mergePoints.H"
#include "fileName.H"
#include "unitConversion.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::typeName
    = "CGALSurfaceMesh";


const std::map<Foam::word, Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::validShape_>
Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::validShapes_ =
{
//    {"plane",      validShape_::PLANE},
    {"box",        validShape_::BOX},
//    {"rotatedBox", validShape_::ROTATEDBOX},
    {"sphere",     validShape_::SPHERE},
    {"cylinder",   validShape_::CYLINDER},
    {"cone",       validShape_::CONE},
//    {"torus",      validShape_::TORUS},
    {"file",       validShape_::FILE}
};


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::createBox
(
    const dictionary& shape
)
{
    const point minBound(shape.get<point>("min"));
    const point maxBound(shape.get<point>("max"));

    surface_.reserve(8, 24, 12);

    surface_.add_vertex(CPoint(minBound.x(), minBound.y(), minBound.z()));
    surface_.add_vertex(CPoint(maxBound.x(), minBound.y(), minBound.z()));
    surface_.add_vertex(CPoint(maxBound.x(), maxBound.y(), minBound.z()));
    surface_.add_vertex(CPoint(minBound.x(), maxBound.y(), minBound.z()));
    surface_.add_vertex(CPoint(minBound.x(), minBound.y(), maxBound.z()));
    surface_.add_vertex(CPoint(maxBound.x(), minBound.y(), maxBound.z()));
    surface_.add_vertex(CPoint(maxBound.x(), maxBound.y(), maxBound.z()));
    surface_.add_vertex(CPoint(minBound.x(), maxBound.y(), maxBound.z()));

    surface_.add_face(CMVIndex(0), CMVIndex(3), CMVIndex(1));
    surface_.add_face(CMVIndex(1), CMVIndex(3), CMVIndex(2));
    surface_.add_face(CMVIndex(1), CMVIndex(6), CMVIndex(5));
    surface_.add_face(CMVIndex(1), CMVIndex(2), CMVIndex(6));
    surface_.add_face(CMVIndex(0), CMVIndex(4), CMVIndex(3));
    surface_.add_face(CMVIndex(3), CMVIndex(4), CMVIndex(7));
    surface_.add_face(CMVIndex(0), CMVIndex(1), CMVIndex(4));
    surface_.add_face(CMVIndex(1), CMVIndex(5), CMVIndex(4));
    surface_.add_face(CMVIndex(2), CMVIndex(3), CMVIndex(6));
    surface_.add_face(CMVIndex(3), CMVIndex(7), CMVIndex(6));
    surface_.add_face(CMVIndex(4), CMVIndex(5), CMVIndex(6));
    surface_.add_face(CMVIndex(4), CMVIndex(6), CMVIndex(7));
}


void Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::createSphere
(
    const dictionary& shape
)
{
    const scalar radius(shape.get<scalar>("radius"));
    const point centre(shape.get<point>("centre"));
    const label subdivision(shape.get<label>("subdivision"));

    // Initial points and faces
    //   from PyMesh (https://github.com/PyMesh/PyMesh)
    const scalar paraR(0.5 * (1.0 + Foam::sqrt(5.0)));
    pointField paraPoints
    (
        std::initializer_list<point>
        {
            point(  -1.0,  paraR,    0.0), point(   1.0,  paraR,    0.0),
            point(  -1.0, -paraR,    0.0), point(   1.0, -paraR,    0.0),
            point(   0.0,   -1.0,  paraR), point(   0.0,    1.0,  paraR),
            point(   0.0,   -1.0, -paraR), point(   0.0,    1.0, -paraR),
            point( paraR,    0.0,   -1.0), point( paraR,    0.0,    1.0),
            point(-paraR,    0.0,   -1.0), point(-paraR,    0.0,    1.0)
        }
    );

    faceList paraFaces
    (
        std::initializer_list<face>
        {
            face(triFace( 0, 11,  5)), face(triFace( 0,  5,  1)),
            face(triFace( 0,  1,  7)), face(triFace( 0,  7, 10)),
            face(triFace( 0, 10, 11)), face(triFace( 1,  5,  9)),
            face(triFace( 5, 11,  4)), face(triFace(11, 10,  2)),
            face(triFace(10,  7,  6)), face(triFace( 7,  1,  8)),
            face(triFace( 3,  9,  4)), face(triFace( 3,  4,  2)),
            face(triFace( 3,  2,  6)), face(triFace( 3,  6,  8)),
            face(triFace( 3,  8,  9)), face(triFace( 5,  4,  9)),
            face(triFace( 2,  4, 11)), face(triFace( 6,  2, 10)),
            face(triFace( 8,  6,  7)), face(triFace( 9,  8,  1))
        }
    );

    // Subdivision
    for (label subdivLevel = 1; subdivLevel <= subdivision; subdivLevel++)
    {
        // Copy previous faces and clear
        //const pointField paraPointsCopy(paraPoints);
        const faceList paraFacesCopy(paraFaces);
        //paraPoints.clear();
        paraFaces.clear();

        // Perform subdivision for each triangle
        forAll (paraFacesCopy, facei)
        {
            const face& fa(paraFacesCopy[facei]);

            /* Compute 3 extra points by spliting half on each edge
            //           P1
            //          / \
            //   newP1 *---* newP3
            //        / \ / \
            //      P2---*---P3
            //         newP2
            */
            const point& P1(paraPoints[fa[0]]);
            const point& P2(paraPoints[fa[1]]);
            const point& P3(paraPoints[fa[2]]);

            const point newP1(0.5 * (P1 + P2));
            const point newP2(0.5 * (P3 + P2));
            const point newP3(0.5 * (P1 + P3));

            const label maxPtId(paraPoints.size() - 1);

            paraPoints.append(newP1);
            paraPoints.append(newP2);
            paraPoints.append(newP3);

            paraFaces.append(face(triFace(    fa[0], maxPtId+1, maxPtId+3)));
            paraFaces.append(face(triFace(maxPtId+1,     fa[1], maxPtId+2)));
            paraFaces.append(face(triFace(maxPtId+1, maxPtId+2, maxPtId+3)));
            paraFaces.append(face(triFace(maxPtId+3, maxPtId+2,     fa[2])));
        }
    }

    // Merge duplicated points and update point labels of faces
    labelList pointToUnique;
    inplaceMergePoints<pointField>(paraPoints, 1e-8, false, pointToUnique);

    forAll (paraFaces, facei)
    {
        face& fa(paraFaces[facei]);

        fa[0] = pointToUnique[fa[0]];
        fa[1] = pointToUnique[fa[1]];
        fa[2] = pointToUnique[fa[2]];
    }

    // Scaling
    for (point& pointi: paraPoints)
    {
        pointi = pointi * radius / mag(pointi) + centre;
    }

    // Convert to CGAL surface mesh
    convert(paraPoints, paraFaces);
}


void Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::createCylinder
(
    const dictionary& shape
)
{
    const point point1(shape.get<point>("point1"));
    const point point2(shape.get<point>("point2"));
    const scalar radius(shape.get<scalar>("radius"));
    const label nvertices(shape.get<label>("vertices"));
    const scalar innerRadius(shape.getOrDefault<scalar>("innerRadius", 0.0));

    dictionary createCylinderDict;

    createCylinderDict.add("point1", point1);
    createCylinderDict.add("point2", point2);
    createCylinderDict.add("radius1", radius);
    createCylinderDict.add("radius2", radius);
    createCylinderDict.add("vertices", nvertices);
    createCylinderDict.add("innerRadius1", innerRadius);
    createCylinderDict.add("innerRadius2", innerRadius);

    createCone(createCylinderDict);
}


void Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::createCone
(
    const dictionary& shape
)
{
    const point point1(shape.get<point>("point1"));
    const point point2(shape.get<point>("point2"));
    const scalar radius1(shape.get<scalar>("radius1"));
    const scalar radius2(shape.get<scalar>("radius2"));
    const label nvertices(shape.get<label>("vertices"));
    const scalar innerRadius1(shape.getOrDefault<scalar>("innerRadius1", 0.0));
    const scalar innerRadius2(shape.getOrDefault<scalar>("innerRadius2", 0.0));

    pointField conePoints(label(nvertices*4), point::zero);
    faceList coneFaces(label(nvertices*8), face(triFace(0,0,0)));

    const scalar deltaTheta
    (
        Foam::constant::mathematical::twoPi / scalar(nvertices)
    );

    vector normal(point2 - point1);
    normal /= (mag(normal) + SMALL);

    vector v1(vector(normal.z(), 0.0, -normal.x()));
    v1 /= (mag(v1) + SMALL);

    const vector v2(normal ^ v1);

    conePoints.clear();
    coneFaces.clear();

    for (label i = 0; i < nvertices; i++)
    {
        const scalar thetai(scalar(i) * deltaTheta);
        const scalar cosThetai(Foam::cos(thetai));
        const scalar sinThetai(Foam::sin(thetai));
        const point circlei(cosThetai * v1 + sinThetai * v2);

        conePoints.append(point1 + radius1 * circlei);
        conePoints.append(point1 + innerRadius1 * circlei);
        conePoints.append(point2 + radius2 * circlei);
        conePoints.append(point2 + innerRadius2 * circlei);
    }

    for (label i = 0; i < nvertices; i++)
    {
        const label nexti((i+1)%nvertices);
        
        coneFaces.append(face(triFace(      4*i,   4*nexti,     4*i+2)));
        coneFaces.append(face(triFace(  4*nexti, 4*nexti+2,     4*i+2)));

        coneFaces.append(face(triFace(    4*i+1,     4*i+3, 4*nexti+1)));
        coneFaces.append(face(triFace(4*nexti+1,     4*i+3, 4*nexti+3)));

        coneFaces.append(face(triFace(      4*i,     4*i+1,   4*nexti)));
        coneFaces.append(face(triFace(  4*nexti,     4*i+1, 4*nexti+1)));

        coneFaces.append(face(triFace(    4*i+2, 4*nexti+2,     4*i+3)));
        coneFaces.append(face(triFace(4*nexti+2, 4*nexti+3,     4*i+3)));
    }

    // Merge duplicated points and update point labels of faces
    labelList pointToUnique;
    inplaceMergePoints<pointField>(conePoints, 1e-8, false, pointToUnique);

    forAll (coneFaces, facei)
    {
        face& fa(coneFaces[facei]);

        fa[0] = pointToUnique[fa[0]];
        fa[1] = pointToUnique[fa[1]];
        fa[2] = pointToUnique[fa[2]];
    }

    // Convert to CGAL surface mesh
    convert(conePoints, coneFaces);
}


void Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::createFromFile
(
    const dictionary& shape
)
{
    const fileName shapePath(shape.get<fileName>("path"));

    if (!exists(shapePath))
    {
        Info<< "        Warning: file "
                << shapePath
                << " in "
                << shape.relativeName(true)
                << " does not exist. Skip."
                << nl << endl;

        return;
    }

    meshedSurface mSurf(shapePath);
    mSurf.scalePoints(shape.getOrDefault<scalar>("scale", 1.0));

    // Merge duplicated points and update point labels of faces
    pointField localPoints(mSurf.points());
    faceList localFaces(mSurf.surfFaces());
    labelList pointToUnique;
    inplaceMergePoints<pointField>(localPoints, 1e-8, false, pointToUnique);

    forAll (localFaces, facei)
    {
        face& fa(localFaces[facei]);

        fa[0] = pointToUnique[fa[0]];
        fa[1] = pointToUnique[fa[1]];
        fa[2] = pointToUnique[fa[2]];
    }

    point centerPoint(point::zero);
    forAll(localPoints, pointi)
    {
        centerPoint += localPoints[pointi];
    }
    centerPoint /= scalar(localPoints.size());

    // Move to origin O for rotation
    localPoints -= centerPoint;

    const vector rotateVec(shape.getOrDefault<vector>("rotate", vector::zero));
    const scalar rAlpha(Foam::degToRad(rotateVec.x()));
    const scalar rBeta(Foam::degToRad(rotateVec.y()));
    const scalar rGamma(Foam::degToRad(rotateVec.z()));

    const tensor rotateMatrix
    (
        Foam::cos(rBeta)*Foam::cos(rGamma),
        Foam::sin(rAlpha)*Foam::sin(rBeta)*Foam::cos(rGamma)-Foam::cos(rAlpha)*Foam::sin(rGamma),
        Foam::cos(rAlpha)*Foam::sin(rBeta)*Foam::cos(rGamma)+Foam::sin(rAlpha)*Foam::sin(rGamma),
        Foam::cos(rBeta)*Foam::sin(rGamma),
        Foam::sin(rAlpha)*Foam::sin(rBeta)*Foam::sin(rGamma)+Foam::cos(rAlpha)*Foam::cos(rGamma),
        Foam::cos(rAlpha)*Foam::sin(rBeta)*Foam::sin(rGamma)-Foam::sin(rAlpha)*Foam::cos(rGamma),
        -Foam::sin(rBeta),
        Foam::sin(rAlpha)*Foam::cos(rBeta),
        Foam::cos(rAlpha)*Foam::cos(rBeta)
    );

    localPoints = rotateMatrix & localPoints;

    // Move back to center point
    localPoints += centerPoint;


    localPoints += shape.getOrDefault<vector>("translate", vector::zero);

    // Convert to CGAL surface mesh
    convert(localPoints, localFaces);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::CGALSurfaceMesh()
:
    surface_(CMesh())
{
    surface_.clear();
}


Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::CGALSurfaceMesh
(
    const dictionary& shape
)
{
    surface_.clear();

    const word shapeType(shape.get<word>("type"));

    std::map<word, validShape_>::const_iterator iter
    (
        validShapes_.find(shapeType)
    );

    if (iter != validShapes_.end())
    {
        switch (iter->second)
        {
            //case validShape_::PLANE:
            //    break;

            case validShape_::BOX:
                createBox(shape);
                break;

            //case validShape_::ROTATEDBOX:
            //    break;

            case validShape_::SPHERE:
                createSphere(shape);
                break;

            case validShape_::CYLINDER:
                createCylinder(shape);
                break;

            case validShape_::CONE:
                createCone(shape);
                break;

            //case validShape_::TORUS:
            //    break;

            case validShape_::FILE:
                createFromFile(shape);
                break;
        }
    }
    else
    {
        Info<< "        Warning: type "
            << fileName(shapeType)
            << " in "
            << shape.relativeName(true)
            << " is not supported. Skip."
            << endl;
    }
}


Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::CGALSurfaceMesh
(
    const fvMesh& mesh,
    const label cellI,
    const word triangulationScheme
)
{
    surface_.clear();

    const cellList& cells(mesh.cells());
    const faceList& faces(mesh.faces());
    const pointField& points(mesh.points());
    //const vectorField& cCtrs(mesh.cellCentres());
    const vectorField& fCtrs(mesh.faceCentres());
    const labelList& own(mesh.faceOwner());

    // Localize cell point labels
    const cell& c(cells[cellI]);
    const labelList globalPointLabels(c.labels(faces));

    pointField cPts(20);
    faceList fList(20);
    cPts.clear();
    fList.clear();

    forAll(globalPointLabels, pointi)
    {
        cPts.append(points[globalPointLabels[pointi]]);
    }

    if (triangulationScheme == word("centroid"))
    {
        forAll(c, facei)
        {
            const face& fa(faces[c[facei]]);

            if (fa.size() == 3)
            {
                face localFa(3);
                localFa[0] = globalPointLabels.find(fa[0]);
                localFa[1] = globalPointLabels.find(fa[1]);
                localFa[2] = globalPointLabels.find(fa[2]);

                fList.append
                (
                    cellI == own[c[facei]] ? localFa : localFa.reverseFace()
                );
            }
            else
            {
                cPts.append(fCtrs[c[facei]]);

                forAll(fa, pointi)
                {
                    const label nextPointi((pointi + 1) % fa.size());

                    face localFa(3);
                    //localFa.clear();

                    localFa[0] = cPts.size()-1;
                    localFa[1] = globalPointLabels.find(fa[pointi]);
                    localFa[2] = globalPointLabels.find(fa[nextPointi]);

                    fList.append
                    (
                        cellI == own[c[facei]] ? localFa : localFa.reverseFace()
                    );
                }
            }
        }
    }
    else
    {
        forAll(c, facei)
        {
            const face& fa(faces[c[facei]]);

            for (label i = 0; i < fa.size()-2; i++)
            {
                const label nextPointi(i + 1);
                const label nextNextPointi(i + 2);

                face localFa(3);
                localFa[0] = globalPointLabels.find(fa[0]);
                localFa[1] = globalPointLabels.find(fa[i+1]);
                localFa[2] = globalPointLabels.find(fa[i+2]);

                fList.append
                (
                    cellI == own[c[facei]] ? localFa : localFa.reverseFace()
                );
            }
        }
    }

    // Convert to CGAL surface mesh
    convert(cPts, fList);
}


Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::CGALSurfaceMesh
(
    const DynamicList<point>& dymPointList,
    const DynamicList<face>& dymfaceList
)
{
    pointField pField(dymPointList.size() + dymfaceList.size());
    pField.clear();
    forAll(dymPointList, pointI)
    {
        pField.append(dymPointList[pointI]);
    }

    DynamicList<face> fList(10*dymfaceList.size());
    fList.clear();
    forAll(dymfaceList, faceI)
    {
        const face& fa(dymfaceList[faceI]);

        if (fa.size() > 3)
        {
            const point fCentre(fa.centre(dymPointList));
            pField.append(fCentre);

            forAll(fa, pointI)
            {
                fList.append
                (
                    face
                    (
                        triFace
                        (
                            pField.size()-1,
                            fa.thisLabel(pointI),
                            fa.nextLabel(pointI)
                        )
                    )
                );
            }
        }
        else
        {
            fList.append(fa);
        }
    }

    convert(pField, fList);
}


Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::CGALSurfaceMesh
(
    const meshedSurface& mSurf
)
{
    surface_.clear();

    convert(mSurf.points(), mSurf.surfFaces());
}


Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::CGALSurfaceMesh
(
    const CMesh& cgalSurfMesh
)
:
    surface_(cgalSurfMesh)
{
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::scalar Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::volume() const
{
    if (surface_.number_of_vertices() > 0)
    {
        return scalar(PMP::volume(surface_));
    }
    else
    {
        return 0.0;
    }
}


void Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::convert
(
    const pointField& cPts,
    const faceList& fList
)
{
    surface_.clear();

    unsigned int nV(cPts.size());
    unsigned int nF(fList.size());
    unsigned int nE(cPts.size() + fList.size() - 2);

    surface_.reserve(nV, (std::max)(3*nV, nE), nF);
    //surface_.reserve(nV, 3*nF, nF);

    forAll(cPts, ptI)
    {
        const point& pt(cPts[ptI]);
        surface_.add_vertex(CPoint(pt.x(), pt.y(), pt.z()));
    }

    forAll(fList, faI)
    {
        const face& fa(fList[faI]);
        surface_.add_face(CMVIndex(fa[0]), CMVIndex(fa[1]), CMVIndex(fa[2]));
    }
}


Foam::pointField Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::getPointField() const
{
    const label nPoints(surface_.number_of_vertices());

    pointField points;

    for (label ptI = 0; ptI < nPoints; ptI++)
    {
        const CPoint& CPt(surface_.point(CMVIndex(ptI)));

        points.append
        (
            point(scalar(CPt.x()), scalar(CPt.y()), scalar(CPt.z()))
        );
    }

    return points;
}


Foam::faceList Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::getFaceList() const
{
    faceList faces;

    for
    (
        CMFIter fIt = surface_.faces_begin();
        fIt != surface_.faces_end();
        ++fIt
    )
    {
        face fa;

        CMVFCirt fVIt(surface_.halfedge(*fIt), surface_);
        CMVFCirt fVEnd(fVIt);

        do
        {
            CMSzType CPtI(*fVIt);
            fa.append(label(CPtI));
        }
        while (++fVIt != fVEnd);

        faces.append(fa);
    }

    return faces;
}


bool Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::clip
(
    const vector& normal,
    const scalar distance
)
{
    //PMP::remove_degenerate_faces(surface_);

    return PMP::clip
    (
        surface_,
        CPlane(normal.x(), normal.y(), normal.z(), distance),
        PMP::parameters::clip_volume(true).use_compact_clipper(true).throw_on_self_intersection(true)
    );
}


void Foam::geometricVofExt::CGALVof::CGALSurfaceMesh::write(const fileName& fName) const
{
    const meshedSurface mSurf(getPointField(), getFaceList());

    mSurf.write(fName);
}


// ************************************************************************* //