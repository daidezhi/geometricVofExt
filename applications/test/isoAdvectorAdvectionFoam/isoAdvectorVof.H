{
    // Create VOF solver
    isoAdvection advector(alpha1, phi, U);

    // Reconstruct and write interface at the initial time (t = 0)
    runTime.functionObjects().execute();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "alphaSuSp.H"
        advector.advect(Sp, Su);   //Interface advection

        #include "updateU.H"

        Info<< nl << "Phase-1 volume = "
            << gSum(alpha1.primitiveField() * mesh.V())
            << nl << endl;

        Info<< "Times consumed by steps:" << nl
            << "  Reconstruction  = " << advector.reconstructionTime() << nl
            << "  Advection       = " << advector.advectionTime() << nl
            << endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }
}
