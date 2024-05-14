IOobject dictIO
(
    IOobject::selectIO
    (
        IOobject
        (
            "setVofFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ),
        args.getOrDefault<fileName>("dict", "")
    )
);

IOdictionary setVofFieldDict(dictIO);

Info<< "\nReading " << setVofFieldDict.name() << nl << endl;


const word fieldName(setVofFieldDict.get<word>("field"));

const bool extraBooleanCheck
(
    setVofFieldDict.getOrDefault<bool>("extraBooleanCheck", false)
);

const scalar initialValue
(
    setVofFieldDict.getOrDefault<scalar>("initial", 0.0)
);

const word triangulationScheme
(
    setVofFieldDict.getOrDefault<word>
    (
        "triangulationScheme", "fan"
    )
);

const dictionary& shapes(setVofFieldDict.subDict("shapes"));
