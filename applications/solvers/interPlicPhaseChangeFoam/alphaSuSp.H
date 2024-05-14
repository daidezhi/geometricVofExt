//zeroField Su;
//zeroField Sp;

volScalarField::Internal divU
(
    mesh.moving() ? fvc::div(phi + mesh.phi()) : fvc::div(phi)
);

Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
const volScalarField& vDotcAlphal = vDotAlphal[0]();
const volScalarField& vDotvAlphal = vDotAlphal[1]();
const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

volScalarField::Internal Sp((vDotvmcAlphal+divU)());
volScalarField::Internal Su(vDotcAlphal);