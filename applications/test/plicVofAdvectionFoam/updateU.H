{
    mesh.C();
    mesh.Cf();
    mesh.Sf();

    Info<< "Calculating field U" << endl;

    const scalarField x(mesh.C().component(vector::X));
    const scalarField y(mesh.C().component(vector::Y));
    const scalarField z(mesh.C().component(vector::Z));

    vectorField& Uc = U.primitiveFieldRef();
    Uc.replace(vector::X, 2*sqr(sin(M_PI*x))*sin(2*M_PI*y)*sin(2*M_PI*z));
    Uc.replace(vector::Y, -sin(2*M_PI*x)*sqr(sin(M_PI*y))*sin(2*M_PI*z));
    Uc.replace(vector::Z, -sin(2*M_PI*x)*sin(2*M_PI*y)*sqr(sin(M_PI*z)));

    U.correctBoundaryConditions();


    Info<< "Calculating face flux field phi\n" << endl;

    const vectorField Cf(mesh.Cf().primitiveField());
    const scalarField Xf(Cf.component(vector::X));
    const scalarField Yf(Cf.component(vector::Y));
    const scalarField Zf(Cf.component(vector::Z));
    vectorField Uf(Xf.size());
    Uf.replace(0, 2*sqr(sin(M_PI*Xf))*sin(2*M_PI*Yf)*sin(2*M_PI*Zf));
    Uf.replace(1, -sin(2*M_PI*Xf)*sqr(sin(M_PI*Yf))*sin(2*M_PI*Zf));
    Uf.replace(2, -sin(2*M_PI*Xf)*sin(2*M_PI*Yf)*sqr(sin(M_PI*Zf)));

    scalarField& phic = phi.primitiveFieldRef();
    const vectorField& Sfc = mesh.Sf().primitiveField();
    phic = Uf & Sfc;

    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();
    const surfaceVectorField::Boundary& Sfbf =
        mesh.Sf().boundaryField();
    const surfaceVectorField::Boundary& Cfbf =
        mesh.Cf().boundaryField();

    forAll(phibf, patchi)
    {
        scalarField& phif = phibf[patchi];
        const vectorField& Sff = Sfbf[patchi];
        const vectorField& Cff = Cfbf[patchi];
        const scalarField xf(Cff.component(vector::X));
        const scalarField yf(Cff.component(vector::Y));
        const scalarField zf(Cff.component(vector::Z));
        vectorField Uf(xf.size());
        Uf.replace(0, 2*sqr(sin(M_PI*xf))*sin(2*M_PI*yf)*sin(2*M_PI*zf));
        Uf.replace(1, -sin(2*M_PI*xf)*sqr(sin(M_PI*yf))*sin(2*M_PI*zf));
        Uf.replace(2, -sin(2*M_PI*xf)*sin(2*M_PI*yf)*sqr(sin(M_PI*zf)));
        phif = Uf & Sff;
    }

    const scalar t(runTime.time().value());
    const scalar dt(runTime.deltaT().value());

    if ( period > 0.0 )
    {
        const scalar uFactor
        (
            0.5*(Foam::cos(2.0*M_PI*t/period)
               + Foam::cos(2.0*M_PI*(t + dt)/period))
        );

        phi *= uFactor;
        U *= uFactor;
    }
}
