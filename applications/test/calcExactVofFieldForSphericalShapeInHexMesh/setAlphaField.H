{
    Info<< "Set field " << fieldName << nl << endl;

    const scalarField& cellVolumes(mesh.cellVolumes());

    forAll(submergedCells, cellI)
    {
        alpha1In[submergedCells[cellI]] = scalar(1.0);
    }


    #pragma omp parallel for
    forAll(intersectedCells, cellI)
    {
        const label cellID(intersectedCells[cellI]);
        
        alpha1In[cellID] =
            cutVolume(mesh, cellID, shape) / cellVolumes[cellID];
    }


    alpha1.correctBoundaryConditions();

    Info<< "sum(alpha*V): "
        << gSum(mesh.V() * alpha1In)
        << nl
        << endl;
}
