Foam::word regionName(Foam::polyMesh::defaultRegion);

{
    Foam::Info << "Create mesh";
    if (args.readIfPresent("region", regionName))
    {
        Foam::Info << ' ' << regionName;
    }
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;
}

autoPtr<dynamicFvMesh> meshPtr
(
    dynamicFvMesh::New
    (
        Foam::IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    )
);

dynamicFvMesh& mesh(meshPtr());

mesh.init(true);   // Initialise all (lower levels and current)

Foam::Info << Foam::endl;

//