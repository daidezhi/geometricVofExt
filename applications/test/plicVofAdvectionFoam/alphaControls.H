const dictionary& alphaControls(mesh.solverDict(alpha1.name()));

//label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

const scalar period(alphaControls.lookupOrDefault<scalar>("period", 0.0));
//scalar reverseTime = alphaControls.lookupOrDefault<scalar>("reverseTime", 0.0);
