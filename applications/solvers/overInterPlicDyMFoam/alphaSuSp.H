zeroField Su;
//zeroField Sp;

volScalarField::Internal divU
(
    mesh.moving() ? fvc::div(phi + mesh.phi()) : fvc::div(phi)
);

divU *= interpolatedCells*cellMask;

volScalarField::Internal Sp(divU);
