    // Disable and add some options
    argList::noParallel();          // Disable parallel function
    argList::noFunctionObjects();   // Don't use function objects

    argList::addOption
    (
        "centre",
        "point",
        "Centre of the spherical shape (Default: (0 0 0))"
    );

    argList::addOption
    (
        "radius",
        "scalar",
        "Radius of the spherical shape (Default: 1)"
    );

    argList::addOption
    (
        "field",
        "word",
        "Name of alpha field (Default: alpha.water.exact)"
    );

    argList::addOption
    (
        "targetField",
        "word",
        "Target field for shape error estimation"
    );

    argList::addOption("np", "number", "Number of OpenMP threads");
