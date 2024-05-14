// Prepare file writing
const fileName postDir(args.globalPath() / functionObject::outputPrefix);

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
    << string("E_v(t)").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << string("min(alpha)").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << string("1-max(alpha)").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << string("E_p(t)").c_str();
os  << Foam::setf(ios_base::right) << Foam::setw(numberWidth)
    << string("E_s(t)").c_str()
    << nl;

os.flush();
