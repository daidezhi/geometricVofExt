// Errors
scalar volumeConservationErr(0.0);
scalar minAlpha(0.0);
scalar oneMinusMaxAlpha(0.0);
scalar sharpnessErr(0.0);
scalar shapeErr(0.0);


// Volume conservation
volumeConservationErr = (gSum(alpha1In * cellVol) - exactInitialVol) / exactInitialVol;

// Boundedness
minAlpha = gMin(alpha1In);
oneMinusMaxAlpha = scalar(1.0) - gMax(alpha1In);

// Sharpness
scalarList isoValues(std::initializer_list<scalar>({0.01, 0.99}));

isoSurfaceVolume calcIsoVolumes(alpha1, isoValues);
const scalarList& isoVolumes(calcIsoVolumes.volumes());

isoSurfaceVolume calcIsoVolumesExact(alpha1Ext, isoValues);
const scalarList& isoVolumesExact(calcIsoVolumesExact.volumes());

scalar deltaVol(isoVolumes[0] - isoVolumes[1]);
scalar deltaVolExact(isoVolumesExact[0] - isoVolumesExact[1]);

sharpnessErr = mag(deltaVol - deltaVolExact) / (deltaVolExact + SMALL);

// Shape preservation
shapeErr = gSum(mag(alpha1In - alpha1ExtIn) * cellVol) / gSum(alpha1ExtIn * cellVol);


// Write to file
os  << Foam::setw(2) << "  ";
os  << Foam::setf(ios_base::right) << Foam::setw(5)
    << runTime.timeName();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << volumeConservationErr;
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << minAlpha;
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << oneMinusMaxAlpha;
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << sharpnessErr;
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << shapeErr
    << nl;

os.flush();


Info<< "E_v(t)     = " << volumeConservationErr << endl
    << "alphaMin   = " << minAlpha << endl
    << "1-alphaMax = " << oneMinusMaxAlpha << endl
    << "E_p(t)     = " << sharpnessErr << endl
    << "E_s(t)     = " << shapeErr << endl
    << nl << endl;
