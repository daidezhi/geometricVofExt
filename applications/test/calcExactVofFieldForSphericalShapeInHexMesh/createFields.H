    Info<< "Reading field " << fieldName << nl << endl;

    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    scalarField& alpha1In(alpha1.ref());
    alpha1In = scalar(0.0);
    alpha1.correctBoundaryConditions();


    Info<< "Creating the spherical shape " << nl << endl;

    Sphere shape{vector_t(centre.x(), centre.y(), centre.z()), radius};


    //- Point locations
    //      -1: bounded (submerged)
    //       0: on boundary (on the surface mesh boundary)
    //       1: unbounded (not submerged)
    //      -2: not determined
    labelList pointLocations(mesh.nPoints(), -2);


    //- Submerged cell list
    DynamicList<label> submergedCells(label(0.4*mesh.nCells()));

    //- Intersected cell list
    DynamicList<label> intersectedCells(label(0.4*mesh.nCells()));
