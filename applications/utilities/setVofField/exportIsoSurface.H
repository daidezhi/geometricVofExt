if (setVofFieldDict.getOrDefault<bool>("exportIsoSurface", false))
{
    Info<< "Exporting ISO interface (ISO = 0.5) ..." << endl;

    volPointInterpolation vpi(mesh);

    pointScalarField alpha1P(vpi.interpolate(alpha1));

    if (!exists(outputDir)) {mkDir(outputDir);}

    isoSurfaceTopo isoSurface(mesh, alpha1, alpha1P, 0.5);

    isoSurface.write(fileName(outputDir / "isoSurface.vtp"));

    Info << endl;
}
