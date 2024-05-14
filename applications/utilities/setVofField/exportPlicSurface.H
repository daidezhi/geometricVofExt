if (setVofFieldDict.getOrDefault<bool>("exportPlicSurface", false))
{
    Info<< "Exporting PLIC interface ..." << endl;

    dictionary reconstructionDict;

    reconstructionDict.add
    (
        "mixedCellTol",
        setVofFieldDict.getOrDefault<scalar>("snapTol", 1e-8)
    );

    reconstructionDict.add("orientationMethod", word("RDF"));
    reconstructionDict.add("tol", 1e-6);
    reconstructionDict.add("relTol", 0.01);
    reconstructionDict.add("iterations", 100);

    reconstructionDict.add("mapAlphaField", true);

    geometricVofExt::SimPLIC::reconstruction
        reconstructor(alpha1, reconstructionDict);

    reconstructor.reconstruct();

    if (!exists(outputDir)) {mkDir(outputDir);}

    reconstructor.interface().write(fileName(outputDir / "plicSurface.vtp"));

    Info << endl;
}
