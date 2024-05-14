IOobject dictIO
(
    IOobject::selectIO
    (
        IOobject
        (
            "setShapedFieldsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ),
        args.getOrDefault<fileName>("dict", "")
    )
);

IOdictionary setShapedFieldsDict(dictIO);

Info<< "\nReading " << setShapedFieldsDict.name() << nl << endl;

const bool extraBooleanCheck
(
    setShapedFieldsDict.getOrDefault<bool>("extraBooleanCheck", false)
);

const word triangulationScheme
(
    setShapedFieldsDict.getOrDefault<word>
    (
        "triangulationScheme", "fan"
    )
);


const dictionary& shapes(setShapedFieldsDict.subDict("shapes"));
