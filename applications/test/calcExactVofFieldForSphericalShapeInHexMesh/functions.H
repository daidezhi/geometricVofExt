scalar cutVolume
(
    const dynamicFvMesh& mesh,
    const label cellI,
    const Sphere& shape
)
{
    const pointField& points(mesh.points());
    const labelList& ptList(mesh.cellShapes()[cellI]);

    List<vector_t> ptCoordList(8);

    forAll(ptList, pti)
    {
        const point& pt(points[ptList[pti]]);

        ptCoordList[pti] = vector_t(pt.x(), pt.y(), pt.z());
    }

    Hexahedron hexI
    {
        ptCoordList[0], ptCoordList[1], ptCoordList[2], ptCoordList[3],
        ptCoordList[4], ptCoordList[5], ptCoordList[6], ptCoordList[7]
    };

    return scalar(overlap(shape, hexI));
}
