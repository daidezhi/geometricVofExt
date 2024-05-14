// Reference shape path and dict
const dictionary& alphaSolveDict(mesh.solverDict(alpha1.name()));
const fileName refShapePath(alphaSolveDict.get<fileName>("refShapePath"));
dictionary refShapeDict;
refShapeDict.add("type", "file");
refShapeDict.add("path", fileName(refShapePath / runTime.timeName() / "shape.stlb"));
const geometricVofExt::CGALVof::CGALSurfaceMesh refShape(refShapeDict);


// Reconstruction data
const DynamicList<label>& mixedCells(reconstructor.mixedCells());
const vectorField& normal(reconstructor.interfaceN().primitiveField());
const scalarField& distance(reconstructor.interfaceD().primitiveField());
geometricVofExt::SimPLIC::cutCell& cellCutter(reconstructor.cellCutter());


// Alpha reference
const scalarField& vol(mesh.cellVolumes());
const scalarField& alphaIn(alpha1.primitiveField());


List<geometricVofExt::CGALVof::CGALSurfaceMesh> sSubCellMeshes(mixedCells.size());
List<scalar> subCellVolumes(mixedCells.size());

forAll (mixedCells, i)
{
    // Get subcell geometry
    const label celli(mixedCells[i]);
    const vector& normali(normal[celli]);
    const scalar& distancei(distance[celli]);
    const label cellStatusi(cellCutter.calcSubCell(celli, normali, distancei, false));

    sSubCellMeshes[i]
    =
    geometricVofExt::CGALVof::CGALSurfaceMesh(cellCutter.subCellPoints(), cellCutter.subCellFaces());

    subCellVolumes[i] = cellCutter.volumeOfSubCell();
}


#pragma omp parallel for
forAll (mixedCells, i)
{
    const label celli(mixedCells[i]);

    const geometricVofExt::CGALVof::CGALSurfaceMesh& sSubCellMeshI(sSubCellMeshes[i]);

    geometricVofExt::CGALVof::surfaceMeshBooleanOperation cutSubCellIByShapeMesh
    (
        sSubCellMeshI,
        refShape,
        word("intersection")
    );
    cutSubCellIByShapeMesh.compBooleanOp();

    if (cutSubCellIByShapeMesh.isValid())
    {
        const scalar commonVol(cutSubCellIByShapeMesh.volume());
        symmDiffVol[celli] = 2.0 * (subCellVolumes[i] - commonVol);
    }
}
