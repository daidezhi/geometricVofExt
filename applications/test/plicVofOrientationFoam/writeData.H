// Prepare file writing
const fileName postDir(args.globalPath() / functionObject::outputPrefix);
mkDir(postDir);

const word fName("performanceData");

autoPtr<OFstream> osPtr;
osPtr.reset(new OFstream(postDir / (fName + ".dat")));
if (!osPtr->good())
{
    FatalIOErrorInFunction(osPtr())
        << "Cannot open file"
        << exit(FatalIOError);
}

Ostream& os(osPtr());

os.setf(ios_base::scientific, ios_base::floatfield);
os.precision(9);
os.width(6+21);

const string head
(
    word(fName) + word("  <-  ") + word(Foam::clock::dateTime())
);

const label numberWidth(18);

os  << Foam::setw(1) << "#";
os  << Foam::setw(1) << ' '
    << Foam::setf(ios_base::left) << head.c_str();
os  << nl;

os  << Foam::setw(1) << "#";
os  << Foam::setw(1) << ' ';
os  << Foam::setf(ios_base::right) << Foam::setw(5)
    << string("Time").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << string("phase1Vol").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << string("symmDiffVol").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth-5)
    << string("nMixedCells").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth-5)
    << string("nCells").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << string("reconsTime").c_str()
    << nl;

os  << Foam::setw(2) << "  ";
os  << Foam::setf(ios_base::right) << Foam::setw(5)
    << runTime.timeName();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << gSum(alphaIn*vol);
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << gSum(symmDiffVol);
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth-5)
    << mixedCells.size();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth-5)
    << mesh.nCells();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << reconstructor.reconstructionTime() << nl;

os.flush();

symmDiffVol.write();