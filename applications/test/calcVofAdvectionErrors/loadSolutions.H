Info<< "Load simulated field alpha1" << endl;
volScalarField alpha1
(
    IOobject
    (
        "alpha.water",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Load exact field alpha1Ext\n" << endl;
volScalarField alpha1Ext
(
    IOobject
    (
        "alpha.water.exact",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


const scalarField& alpha1In(alpha1.internalField());

const scalarField& alpha1ExtIn(alpha1Ext.internalField());

const scalarField& cellVol(mesh.cellVolumes());