{
    const pointField& points(mesh.points());
    const cellList& cells(mesh.cells());
    const cellShapeList& shapes(mesh.cellShapes());

    const scalar TSMALL(10.0*SMALL);


    #pragma omp parallel for
    forAll(points, pointI)
    {
        const point& ptI(points[pointI]);

        const scalar distI(mag(ptI - centre));

        if (distI < TSMALL)
        {
            pointLocations[pointI] = 0;
        }
        else
        {
            pointLocations[pointI] = distI < radius ? -1 : 1;
        }
    }


    submergedCells.clear();
    intersectedCells.clear();

    forAll(cells, cellI)
    {
        const labelList& cellIPoints(shapes[cellI]);

        label nPtsInsideShape(0), nPtsOutsideShape(0), nPtsOnShape(0);

        for (const label pointI : cellIPoints)
        {
            if (pointLocations[pointI] == -1)
            {
                nPtsInsideShape++;
            }
            else if (pointLocations[pointI] == 1)
            {
                nPtsOutsideShape++;
            }
            else
            {
                nPtsOnShape++;
            }
        }

        if ((nPtsInsideShape+nPtsOnShape) == 8)
        {
            submergedCells.append(cellI);
            continue;
        }

        if ((nPtsOutsideShape+nPtsOnShape) == 8)
        {
            continue;
        }

        intersectedCells.append(cellI);
    }
}
