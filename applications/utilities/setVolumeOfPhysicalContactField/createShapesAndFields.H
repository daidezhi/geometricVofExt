    // Create contact cylinders
    Info<< "Creating contact cylinders ..." << endl;

    DynamicList<word> shapeNames;
    DynamicList<word> shapeTypes;
    DynamicList<geometricVofExt::CGALVof::CGALSurfaceMesh> CGALShapes;
    DynamicList<word> shapeModes;

    for (const entry& shapei: contacts.contacts())
    {
        const word shapeIName(shapei.dict().dictName());
        const word shapeIType(shapei.dict().get<word>("type"));
        const word shapeIMode(shapei.dict().getOrDefault<word>("mode", "add"));

        Info<< "    shape: " << shapeIName
            << ", type: " << shapeIType
            << ", mode: " << shapeIMode
            << endl;
        
        std::map<word, validMode>::const_iterator iter
        (
            validModes.find(shapeIMode)
        );

        if (iter == validModes.end())
        {
            Info<< "        Warning: mode "
                << fileName(shapeIMode)
                << " in "
                << shapei.dict().relativeName(true)
                << " is not supported. Skip."
                << nl << endl;

            continue;
        }

        CGALShapes.append
        (
            geometricVofExt::CGALVof::CGALSurfaceMesh(shapei.dict())
        );

        shapeNames.append(shapeIName);

        shapeTypes.append(shapeIType);

        shapeModes.append(shapeIMode);
    }


    // Create Fields
    Info<< "\nCreating zero flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVol/dimTime, 0.0)
    );

    Info<< "Creating zero velocity field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, vector::zero)
    );

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

    alpha1In = initialValue;
    alpha1.correctBoundaryConditions();
