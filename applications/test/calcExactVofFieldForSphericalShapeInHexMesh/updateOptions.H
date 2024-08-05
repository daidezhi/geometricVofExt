    const point centre
    (
        args.found("centre") ? args.get<point>("centre") : point::zero
    );

    const scalar radius
    (
        args.found("radius") ? args.get<scalar>("radius") : 1.0
    );

    const word fieldName
    (
        args.found("field") ?
            args.get<word>("field") : word("alpha.water.exact")
    );
