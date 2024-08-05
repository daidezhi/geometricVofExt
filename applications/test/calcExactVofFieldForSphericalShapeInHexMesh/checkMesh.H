{
    const cellShapeList& shapes(mesh.cellShapes());
    const cellList& cells(mesh.cells());

    forAll(cells, cellI)
    {
        const cellShape& shapeI(shapes[cellI]);

        const label indexI(shapeI.model().index());

        if (cellModel::ref(cellModel::HEX).index() != indexI)
        {
            FatalError
                << "The mesh contains non-hexahedral cells, "
                << "which are not supported by the current utility."
                << nl
                << exit(FatalError);
        }
    }
}
