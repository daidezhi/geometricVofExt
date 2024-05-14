const scalar snapAlphaTol
(
    setVofFieldDict.getOrDefault<scalar>("snapTol", 0.0)
);

if (snapAlphaTol > 0.0)
{
    alpha1 = alpha1
            * pos0(alpha1 - snapAlphaTol)
            * neg0(alpha1 - (1.0 - snapAlphaTol))
            + pos0(alpha1 - (1.0 - snapAlphaTol));
}

alpha1 = min(scalar(1), max(scalar(0), alpha1));

alpha1.correctBoundaryConditions();
