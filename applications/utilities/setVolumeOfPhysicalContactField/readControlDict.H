IOobject dictIO
(
    IOobject::selectIO
    (
        IOobject
        (
            "setVolumeOfPhysicalContactFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        ),
        args.getOrDefault<fileName>("dict", "")
    )
);

IOdictionary setVolumeOfPhysicalContactFieldDict(dictIO);

Info<< "Reading " << setVolumeOfPhysicalContactFieldDict.name() << nl << endl;


const word fieldName(setVolumeOfPhysicalContactFieldDict.get<word>("field"));

const bool extraBooleanCheck
(
    setVolumeOfPhysicalContactFieldDict.getOrDefault<bool>("extraBooleanCheck", false)
);

const scalar initialValue
(
    setVolumeOfPhysicalContactFieldDict.getOrDefault<scalar>("initial", 0.0)
);

const word triangulationScheme
(
    setVolumeOfPhysicalContactFieldDict.getOrDefault<word>
    (
        "triangulationScheme", "fan"
    )
);


// Create contacts
geometricVofExt::contactTable contacts(setVolumeOfPhysicalContactFieldDict);

// Info << contacts.contacts() << endl;