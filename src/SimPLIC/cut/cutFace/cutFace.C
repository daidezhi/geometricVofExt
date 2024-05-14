/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 DHI
    Copyright (C) 2018-2019 Johan Roenby
    Copyright (C) 2019-2020 DLR
    Copyright (C) 2024 Dezhi Dai, Argonne National Laboratory (ANL)
-------------------------------------------------------------------------------
License
    This file is part of geometricVofExt, which is a geometric VOF extension
    to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cutFace.H"
#include "triFace.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::geometricVofExt::SimPLIC::cutFace::calcSubFaceCentreAndArea()
{
    const label nPoints(subFacePoints_.size());

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        subFaceCentre_ = (1.0 / 3.0) * 
            (subFacePoints_[0] + subFacePoints_[1] + subFacePoints_[2]);

        subFaceArea_ = 0.5 * ((subFacePoints_[1] - subFacePoints_[0]) ^
            (subFacePoints_[2] - subFacePoints_[0]));
    }
    else
    {
        vector sumN(vector::zero);
        scalar sumA(0.0);
        vector sumAc(vector::zero);

        // initial guess of centre as average of subFacePoints_
        point fCentre(subFacePoints_[0]);
        for (label pi = 1; pi < nPoints; pi++)
        {
            fCentre += subFacePoints_[pi];
        }
        fCentre /= nPoints;

        // loop sub triangles
        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint(subFacePoints_[(pi + 1) % nPoints]);

            vector c(subFacePoints_[pi] + nextPoint + fCentre);
            vector n
            (
                (nextPoint - subFacePoints_[pi]) ^
                (fCentre - subFacePoints_[pi])
            );
            scalar a(mag(n));

            sumN += n;
            sumA += a;
            sumAc += a * c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA < ROOTVSMALL)
        {
            subFaceCentre_ = fCentre;
            subFaceArea_ = vector::zero;
        }
        else
        {
            subFaceCentre_ = (1.0 / 3.0) * sumAc / sumA;
            subFaceArea_ = 0.5 * sumN;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::geometricVofExt::SimPLIC::cutFace::cutFace
(
    const fvMesh& mesh,
    const scalarField& faceFlatness
)
:
    mesh_(mesh),
    faceFlatness_(faceFlatness),
    subFaceCentre_(point::zero),
    subFaceArea_(vector::zero),
    subFacePoints_(10),
    interfacePoints_(2),
    pointDistances_(10),
    faceStatus_(-1)
{
    clearStorage();
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::label Foam::geometricVofExt::SimPLIC::cutFace::calcSubFace
(
    const label faceI,
    const vector& normal,
    const scalar distance
)
{
    const face& f(mesh_.faces()[faceI]);
    const pointField& points(mesh_.points());

    return calcSubFace(f.points(points), normal, distance);
}


Foam::label Foam::geometricVofExt::SimPLIC::cutFace::calcSubFace
(
    const pointField& fPts,
    const vector& normal,
    const scalar distance
)
{
    clearStorage();

    face f(fPts.size());
    forAll(f, i)
    {
        f[i] = i;
    }

    label nSubmergedPoints(0);
    label firstSubmergedPoint(-1);

    const scalar TSMALL(10.0 * SMALL);

    // Loop face vertices
    forAll(f, i)
    {
        scalar distanceI((fPts[i] & normal) + distance);

        // Lift the vertex slightly if it is very close to the plane
        if (mag(distanceI) < TSMALL)
        {
            distanceI += sign(distanceI) * TSMALL;
        }
        
        pointDistances_.append(distanceI);

        if (distanceI < 0.0)
        {
            nSubmergedPoints++;

            if (firstSubmergedPoint == -1)
            {
                firstSubmergedPoint = i;
            }
        }
    }

    if (nSubmergedPoints == f.size())   // Fully submerged face
    {
        faceStatus_ = -1;

        subFaceCentre_ = f.centre(fPts);
        subFaceArea_ = f.areaNormal(fPts);

        return faceStatus_;
    }
    else if (nSubmergedPoints == 0)     // Fully empty face
    {
        faceStatus_ = 1;

        subFaceCentre_ = point::zero;
        subFaceArea_ = vector::zero;

        return faceStatus_;
    }
    else    // Cut
    {
        faceStatus_ = 0;

        // Loop face and append the cuts
        for
        (
            label i = firstSubmergedPoint;
            i < firstSubmergedPoint + f.size();
            ++i
        )
        {
            const label currentId(i % f.size());
            const label nextId((i + 1) % f.size());

            if (pointDistances_[currentId] < 0) // append submerged vertex
            {
                subFacePoints_.append(fPts[currentId]);
            }

            if      // append intersection point
            (
                (
                    pointDistances_[currentId] * pointDistances_[nextId]
                ) < 0
            )
            {
                const scalar weight
                (
                    pointDistances_[currentId] / 
                    (
                        pointDistances_[currentId] - pointDistances_[nextId]
                    )
                );

                const point cutPoint
                (
                    fPts[currentId] + weight * (fPts[nextId] - fPts[currentId])
                );

                subFacePoints_.append(cutPoint);
                interfacePoints_.append(cutPoint);
            }
        }

        if (subFacePoints_.size() >= 3)
        {
            faceStatus_ = 0;

            calcSubFaceCentreAndArea();
        }
        else
        {
            faceStatus_ = -1;

            subFaceCentre_ = f.centre(fPts);
            subFaceArea_ = f.areaNormal(fPts);
        }

        return faceStatus_;
    }
}


Foam::scalar Foam::geometricVofExt::SimPLIC::cutFace::timeIntegratedFaceFlux
(
    const label faceI,
    const vector& normal,
    const scalar distance,
    const scalar Un0,
    const scalar dt,
    const scalar phi,
    const scalar magSf
)
{
    const scalar TSMALL(10.0 * SMALL);

    if (mag(phi) <= TSMALL)
    {
        return 0.0;
    }

    const face& f(mesh_.faces()[faceI]);
    const pointField& points(mesh_.points());
    const pointField fPts(f.points(points));
    const label nPoints(f.size());
    
    if (mag(Un0 * dt) > TSMALL)     // Interface is not stationary
    {
        // Estimate time of arrival to the face points from their ormal
        // distance to the initial interface and the interface normal velocity
        scalarField pTimes(nPoints);
        forAll(f, i)
        {
            scalar pTimeI(((points[f[i]] & normal) + distance) / Un0);
            pTimes[i] = mag(pTimeI) < TSMALL ? 0.0 : pTimeI;
        }

        scalar dVf(0.0);    // Liquid volume

        if (faceFlatness_[faceI] > (1.0-TSMALL))    // Flat face
        {
            dVf = phi / magSf *
                timeIntegratedArea
                (
                    fPts,
                    normal,
                    distance,
                    pTimes,
                    Un0,
                    dt,
                    magSf
                );
        }
        else    // Warped face, triangular decomposition
        {
            pointField fPtsTri(3);
            const triFace fTri(0,1,2);
            scalarField pTimesTri(3);

            fPtsTri[0] = mesh_.faceCentres()[faceI];

            scalar pTimeTri0(((fPtsTri[0] & normal) + distance) / Un0);
            pTimesTri[0] = mag(pTimeTri0) < TSMALL ? 0.0 : pTimeTri0;

            for (label pi = 0; pi < nPoints; ++pi)
            {
                fPtsTri[1] = fPts[pi];
                pTimesTri[1] = pTimes[pi];

                fPtsTri[2] = fPts[(pi + 1) % nPoints];
                pTimesTri[2] = pTimes[(pi + 1) % nPoints];

                const scalar magSfTri(fTri.mag(fPtsTri));

                const scalar phiTri(phi * magSfTri / magSf);

                dVf += phiTri / magSfTri * 
                    timeIntegratedArea
                    (
                        fPtsTri,
                        normal,
                        distance,
                        pTimesTri,
                        Un0,
                        dt,
                        magSfTri
                    );
            }
        }

        return dVf;
    }
    else    // Un0 is almost zero and interface is treated as stationary
    {
        if (faceFlatness_[faceI] > (1.0-TSMALL))    // Flat face
        {
            calcSubFace(faceI, normal, distance);

            const scalar alphaf(mag(subFaceArea_) / magSf);

            return (phi * dt * alphaf);
        }
        else    // Warped face, triangular decomposition
        {
            pointField fPtsTri(3);
            const triFace fTri(0,1,2);

            fPtsTri[0] = mesh_.faceCentres()[faceI];

            scalar dVf(0.0);    // Liquid volume

            for (label pi = 0; pi < nPoints; ++pi)
            {
                fPtsTri[1] = fPts[pi];
                fPtsTri[2] = fPts[(pi + 1) % nPoints];

                const scalar magSfTri(fTri.mag(fPtsTri));

                const scalar phiTri(phi * magSfTri / magSf);

                calcSubFace(fPtsTri, normal, distance);

                const scalar alphafTri(mag(subFaceArea_) / magSfTri);

                dVf += (phiTri * dt * alphafTri);
            }

            return dVf;
        }
    }
}


Foam::scalar Foam::geometricVofExt::SimPLIC::cutFace::timeIntegratedArea
(
    const pointField& fPts, // face points
    const vector& normal,
    const scalar distance,
    const scalarField& pTimes,
    const scalar Un0,
    const scalar dt,
    const scalar magSf
)
{
    const scalar TSMALL(10.0 * SMALL);

    // Initialize time integrated area returned by this function
    scalar tIntArea(0.0);

    // Finding ordering of vertex points
    const labelList order(Foam::sortedOrder(pTimes));
    const scalar firstTime(pTimes[order.first()]);
    const scalar lastTime(pTimes[order.last()]);

    // Dealing with case where face is not cut by interface during time
    // interval [0, dt] because face was already passed by interface
    if (lastTime <= 0.0)
    {
        // If all face cuttings were in the past and cell is filling up (Un0>0)
        // then face must be full during whole time interval
        tIntArea = magSf * dt * pos0(Un0);
        return tIntArea;
    }

    // Dealing with case where face is not cut by interface during time
    // interval [0, dt] because dt is too small for interface to reach closest
    // face point
    if (firstTime >= dt)
    {
        // If all cuttings are in the future but non of them within [0, dt]
        // then if cell is filling up (Un0 > 0) face must be empty during
        // whole time interval
        tIntArea = magSf * dt * (1.0 - pos0(Un0));
        return tIntArea;
    }

    DynamicList<scalar> sortedTimes(pTimes.size());
    sortedTimes.clear();
    scalar prevTime(0.0);

    scalar subAreaOld(0.0), subAreaNew(0.0), subAreaMid(0.0);

    // Special treatment of first sub time interval
    if (firstTime > 0.0)
    {
        // If firstTime > 0 the face is uncut in the time interval
        // [0, firstTime] and hence fully submerged in fluid A or B.
        // If Un0 > 0 cell is filling up and it must initially be empty.
        // If Un0 < 0 cell must initially be fully immersed in fluid A (liquid).
        subAreaOld = magSf * (scalar(1) - pos0(Un0));
        tIntArea = subAreaOld * firstTime;
        sortedTimes.append(firstTime);
        prevTime = firstTime;
    }
    else
    {
        // If firstTime <= 0 then face is initially cut
        sortedTimes.append(0.0);
        prevTime = 0.0;

        calcSubFace(fPts, normal, distance);
        subAreaOld = mag(subFaceArea_);
    }

    const scalar smallTime(max(TSMALL/mag(Un0), TSMALL));

    forAll(order, ti)
    {
        const scalar timeI(pTimes[order[ti]]);

        if ( timeI > (prevTime + smallTime) && timeI < dt)
        {
            sortedTimes.append(timeI);
            prevTime = timeI;
        }
    }

    if (lastTime > dt)
    {
        sortedTimes.append(dt);
    }
    else
    {
        // Interface will leave the face at lastTime and face will be fully
        // in fluid A or fluid B in the time interval from lastTime to dt.
        tIntArea += magSf * (dt - lastTime) * pos0(Un0);
    }

    for (int k = 0; k < sortedTimes.size()-1; k++)
    {
        const scalar tauOld(sortedTimes[k]);
        const scalar tauNew(sortedTimes[k+1]);
        const scalar deltaTau(0.5 * (tauNew - tauOld));

        calcSubFace(fPts, normal, distance-tauNew*Un0);
        subAreaNew = mag(subFaceArea_);

        calcSubFace(fPts, normal, distance-(tauOld+deltaTau)*Un0);
        subAreaMid = mag(subFaceArea_);

        // Simpson's rule
        tIntArea += (deltaTau/3.0) * (subAreaOld + 4.0*subAreaMid + subAreaNew);

        subAreaOld = subAreaNew;
    }

    return tIntArea;
}


// ************************************************************************* //