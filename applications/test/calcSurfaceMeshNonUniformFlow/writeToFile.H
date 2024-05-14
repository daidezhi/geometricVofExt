// Write to file
autoPtr<OFstream> osPtr;
const word fName(importName.nameLessExt() + "AreasAndVolumes");
const fileName outputDir(args.globalPath() / functionObject::outputPrefix);

osPtr.reset(new OFstream(outputDir / (fName + ".dat")));
if (!osPtr->good())
{
    FatalIOErrorInFunction(osPtr())
        << "Cannot open file"
        << exit(FatalIOError);
}

Ostream& os(osPtr());

os.setf(ios_base::fixed, ios_base::floatfield);
os.precision(15);
os.width(6+21);

const string head
(
    word(fName) + word("  <-  ") + word(Foam::clock::dateTime())
);

os  << Foam::setw(1) << "#";
os  << Foam::setw(1) << ' '
    << Foam::setf(ios_base::left) << head.c_str();
os  << nl;

os  << Foam::setw(1) << "#";
os  << Foam::setw(1) << ' ';
os  << Foam::setf(ios_base::left) << Foam::setw(10)
    << string("Time [s]").c_str();
os  << Foam::setf(ios_base::left) << Foam::setw(20)
    << string("Min AR [-]").c_str();
os  << Foam::setf(ios_base::left) << Foam::setw(25)
    << string("Max AR [-]").c_str();
os  << Foam::setf(ios_base::left) << Foam::setw(25)
    << string("Avg AR [-]").c_str();
os  << Foam::setf(ios_base::left) << Foam::setw(20)
    << string("Area [m^2]").c_str();
os  << Foam::setf(ios_base::left) << Foam::setw(12)
    << string("Volume [m^3]").c_str()
    << nl;

os  << Foam::setw(2) << "  ";
os  << Foam::setf(ios_base::left) << Foam::setw(10) << runTime.timeName();
os  << Foam::setf(ios_base::left) << Foam::setw(20) << gMin(shapeAspectRatios);
os  << Foam::setf(ios_base::left) << Foam::setw(25) << gMax(shapeAspectRatios);
os  << Foam::setf(ios_base::left) << Foam::setw(25) << gAverage(shapeAspectRatios);
os  << Foam::setf(ios_base::left) << Foam::setw(20) << shapeArea;
os  << Foam::setf(ios_base::left) << Foam::setw(17) << shapeVol << nl;


const fileName writePath0
(
    outputDir / importName.nameLessExt() / runTime.timeName()
);

if (!exists(writePath0))
{
    mkDir(writePath0);
}

mSurf0.write(writePath0 / "shape.vtk");
mSurf0.write(writePath0 / "shape.vtp");
mSurf0.write(writePath0 / "shape.stlb");
