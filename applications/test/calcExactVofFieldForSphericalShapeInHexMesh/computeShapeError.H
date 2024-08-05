{
    if (args.found("targetField"))
    {
        const word targetFieldName(args.get<word>("targetField"));

        volScalarField alphaTarget
        (
            IOobject
            (
                targetFieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        const scalarField& alphaExtIn(alpha1.internalField());
        const scalarField& alphaTargetIn(alphaTarget.internalField());
        const scalarField& cellVol(mesh.cellVolumes());

        const scalar shapeErr
        (
            gSum(mag(alphaTargetIn - alphaExtIn) * cellVol)
          / gSum(alphaExtIn * cellVol)
        );

        Info << "Shape error = " << shapeErr << nl << endl;
    }
}
