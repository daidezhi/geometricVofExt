{
    Info<< "Initializing field " << fieldName << nl << endl;

    const scalarField& cellVolumes(mesh.cellVolumes());

    forAll(shapeNames, shapei)
    {
        const Foam::geometricVofExt::CGALVof::CGALSurfaceMesh&
            CGALSMeshI(CGALShapes[shapei]);

        Info<< "    Setting with shape: " << shapeNames[shapei]
            << ", type: " << shapeTypes[shapei]
            << ", volume: " << CGALSMeshI.volume()
            << ", mode: " << shapeModes[shapei]
            << endl;
        
        std::map<word, validMode>::const_iterator iter
        (
            validModes.find(shapeModes[shapei])
        );
        
        Info << "      Calculating mesh point and cell locations ..." << endl;

        const Foam::geometricVofExt::CGALVof::cellToSurfaceMeshLocation cellLocs
        (
            mesh,
            CGALSMeshI,
            extraBooleanCheck
        );
        const DynamicList<label>& intersectedCells(cellLocs.intersectedCells());
        const DynamicList<label>& submergedCells(cellLocs.submergedCells());
        const DynamicList<label>& aridCells(cellLocs.aridCells());

        Info << "      Calculating fraction values ..." << endl;

        const label nIntersectedCells(intersectedCells.size());
        const label nSubmergedCells(submergedCells.size());

        List<scalar> shapeAlpha(nIntersectedCells + nSubmergedCells);

        #pragma omp parallel for
        forAll(intersectedCells, cellI)
        {
            const Foam::geometricVofExt::CGALVof::CGALSurfaceMesh cCell
            (
                mesh,
                intersectedCells[cellI],
                triangulationScheme
            );

            Foam::geometricVofExt::CGALVof::surfaceMeshBooleanOperation booleanOp
            (
                CGALSMeshI,
                cCell,
                word("intersection")
            );

            booleanOp.compBooleanOp();

            shapeAlpha[cellI] =
                booleanOp.volume() / cellVolumes[intersectedCells[cellI]];
        }

        forAll(submergedCells, cellI)
        {
            shapeAlpha[nIntersectedCells + cellI] = 1.0;
        }

        Info << "      Updating fraction values ..." << endl;
        switch (iter->second)
        {
            case validMode::ADD:
                forAll(intersectedCells, cellI)
                {
                    alpha1In[intersectedCells[cellI]] += shapeAlpha[cellI];
                }

                forAll(submergedCells, cellI)
                {
                    alpha1In[submergedCells[cellI]] +=
                        shapeAlpha[nIntersectedCells + cellI];
                }

                break;
            case validMode::SUBTRACT:
                forAll(intersectedCells, cellI)
                {
                    alpha1In[intersectedCells[cellI]] -= shapeAlpha[cellI];
                }

                forAll(submergedCells, cellI)
                {
                    alpha1In[submergedCells[cellI]] -=
                        shapeAlpha[nIntersectedCells + cellI];
                }

                break;
        }

        #include "applyBruteForceBounding.H"
    }

    Info<< nl
        << "sum(alpha*V): "
        << gSum(mesh.V() * alpha1In)
        << nl
        << endl;
}
