#ifndef motionODE_H
#define motionODE_H

#include "ODESystem.H"
#include "ODESolver.H"

#include <cmath>

namespace Foam
{

class motionODE
:
    public ODESystem
{
private:

    scalar period_;

public:

    motionODE(const scalar period)
    :
        period_(period)
    {}

    label nEqns() const
    {
        return 3;
    }

    void derivatives
    (
        const scalar t,         // time
        const scalarField& s,   // motion distance (x y z) = (s[0], s[1], s[2])
        scalarField& dsdt       // velocity
    ) const
    {
        const scalar pi(M_PI);
        const scalar twoPi(2.0 * M_PI);

        dsdt[0] = 2.0*sqr(sin(pi*s[0])) * sin(twoPi*s[1])   * sin(twoPi*s[2]);
        dsdt[1] =      -sin(twoPi*s[0]) * sqr(sin(pi*s[1])) * sin(twoPi*s[2]);
        dsdt[2] =      -sin(twoPi*s[0]) * sin(twoPi*s[1])   * sqr(sin(pi*s[2]));
        dsdt *= cos(twoPi * t / period_);
    }

    void jacobian
    (
        const scalar t,
        const scalarField& s,
        scalarField& dfdt,
        scalarSquareMatrix& dfds
    ) const
    {
        const scalar pi(M_PI);
        const scalar twoPi(2.0 * M_PI);
        const scalar cosTwoPiTime(cos(twoPi * t / period_));

        dfdt[0] = 2.0*sqr(sin(pi*s[0])) * sin(twoPi*s[1])   * sin(twoPi*s[2]);
        dfdt[1] =      -sin(twoPi*s[0]) * sqr(sin(pi*s[1])) * sin(twoPi*s[2]);
        dfdt[2] =      -sin(twoPi*s[0]) * sin(twoPi*s[1])   * sqr(sin(pi*s[2]));
        dfdt *= (-(twoPi / period_) * sin(twoPi * t / period_));

        dfds(0, 0) = 2.0*twoPi*sin(pi*s[0])*cos(pi*s[0]) * sin(twoPi*s[1]) * sin(twoPi*s[2]);
        dfds(0, 1) = 2.0*sqr(sin(pi*s[0])) * twoPi*cos(twoPi*s[1]) * sin(twoPi*s[2]);
        dfds(0, 2) = 2.0*sqr(sin(pi*s[0])) * sin(twoPi*s[1]) * twoPi*cos(twoPi*s[2]);
        dfds(0, 0) *= cosTwoPiTime;
        dfds(0, 1) *= cosTwoPiTime;
        dfds(0, 2) *= cosTwoPiTime;

        dfds(1, 0) = -twoPi*cos(twoPi*s[0]) * sqr(sin(pi*s[1])) * sin(twoPi*s[2]);
        dfds(1, 1) = -sin(twoPi*s[0]) * twoPi*sin(pi*s[1])*cos(pi*s[1]) * sin(twoPi*s[2]);
        dfds(1, 2) = -sin(twoPi*s[0]) * sqr(sin(pi*s[1])) * twoPi*cos(twoPi*s[2]);
        dfds(1, 0) *= cosTwoPiTime;
        dfds(1, 1) *= cosTwoPiTime;
        dfds(1, 2) *= cosTwoPiTime;

        dfds(2, 0) = -twoPi*cos(twoPi*s[0]) * sin(twoPi*s[1]) * sqr(sin(pi*s[2]));
        dfds(2, 1) = -sin(twoPi*s[0]) * twoPi*cos(twoPi*s[1]) * sqr(sin(pi*s[2]));
        dfds(2, 2) = -sin(twoPi*s[0]) * sin(twoPi*s[1]) * twoPi*sin(pi*s[2])*cos(pi*s[2]);
        dfds(2, 0) *= cosTwoPiTime;
        dfds(2, 1) *= cosTwoPiTime;
        dfds(2, 2) *= cosTwoPiTime;
    }
};

}

#endif
