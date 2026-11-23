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

    // Face length scale: largest distance from the first vertex
    scalar magSqrLf(0);
    for (label i = 1; i < fPts.size(); ++i)
    {
        magSqrLf = max(magSqrLf, magSqr(fPts[i] - fPts[0]));
    }

    // Vertex lift relative to the face size (was an absolute 1e-14 m).
    // It never exceeds 1e-14 of the face size, so it cannot distort the
    // geometry. Round-off below it is harmless: the lift keeps the sign,
    // so the vertex classification stays consistent. VSMALL keeps the lift
    // non-zero for degenerate faces.
    const scalar liftTol(TSMALL*Foam::sqrt(magSqrLf) + VSMALL);

    // Loop face vertices
    forAll(f, i)
    {
        scalar distanceI((fPts[i] & normal) + distance);

        // Lift the vertex slightly if it is very close to the plane
        if (mag(distanceI) < liftTol)
        {
            distanceI += sign(distanceI) * liftTol;
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

    // if (mag(phi) < ROOTVSMALL)
    // {
    //     return 0.0;
    // }

    const face& f(mesh_.faces()[faceI]);
    const pointField& points(mesh_.points());
    const pointField fPts(f.points(points));
    const label nPoints(f.size());

    // Arrival-time snapping tolerance (distance), relative to the face size
    const scalar snapTol(TSMALL*Foam::sqrt(magSf) + VSMALL);
    
    // Interface is treated as stationary if it moves less than 1e-14 of the
    // face size during dt (was |Un0| > 1e-12 m/s)
    if (mag(Un0)*dt > TSMALL*Foam::sqrt(magSf))
    {
        // Estimate time of arrival to the face points from their ormal
        // distance to the initial interface and the interface normal velocity
        scalarField pTimes(nPoints);
        forAll(f, i)
        {
            const scalar distI((points[f[i]] & normal) + distance);
            pTimes[i] = mag(distI) < snapTol ? 0.0 : distI/Un0;
        }

        scalar dVf(0.0);    // Liquid volume

        if (faceFlatness_[faceI] > (1.0-TSMALL))    // Flat face
        {
            dVf = phi / magSf *
                timeIntegratedArea
                (
                    faceI,
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

            const scalar dist0((fPtsTri[0] & normal) + distance);
            pTimesTri[0] = mag(dist0) < snapTol ? 0.0 : dist0/Un0;

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
                        faceI,
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
    const label faceI,
    const pointField& fPts, // face points
    const vector& normal,
    const scalar distance,
    const scalarField& pTimes,
    const scalar Un0,
    const scalar dt,
    const scalar magSf
)
{
    //const scalar TSMALL(10.0 * SMALL);

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
        tIntArea = magSf * dt * neg(Un0);
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
        subAreaOld = magSf * neg(Un0);
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

    //const scalar smallTime(max(TSMALL/mag(Un0), TSMALL));
    const scalar smallTime(1e-6*dt);    // was max(1e-6*dt, 1e-14 s)

    forAll(order, ti)
    {
        const scalar timeI(pTimes[order[ti]]);

        if ( timeI > (prevTime + smallTime) && timeI <= dt)
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


// Foam::scalar Foam::geometricVofExt::SimPLIC::cutFace::timeIntegratedArea
// (
//     const label faceI,
//     const pointField& fPts, // face points
//     const vector& normal,
//     const scalar distance,
//     const scalarField& pTimes,
//     const scalar Un0,
//     const scalar dt,
//     const scalar magSf
// )
// {
//     const scalar TSMALL(10.0 * SMALL);

//     // Initialize time integrated area returned by this function
//     scalar tIntArea(0.0);

//     // Finding ordering of vertex points
//     const labelList order(Foam::sortedOrder(pTimes));
//     const scalar firstTime(pTimes[order.first()]);
//     const scalar lastTime(pTimes[order.last()]);

//     // Dealing with case where face is not cut by interface during time
//     // interval [0, dt] because face was already passed by interface
//     if (lastTime <= 0.0)
//     {
//         // If all face cuttings were in the past and cell is filling up (Un0>0)
//         // then face must be full during whole time interval
//         tIntArea = magSf * dt * pos0(Un0);
//         return tIntArea;
//     }

//     // Dealing with case where face is not cut by interface during time
//     // interval [0, dt] because dt is too small for interface to reach closest
//     // face point
//     if (firstTime >= dt)
//     {
//         // If all cuttings are in the future but non of them within [0, dt]
//         // then if cell is filling up (Un0 > 0) face must be empty during
//         // whole time interval
//         tIntArea = magSf * dt * neg(Un0);
//         return tIntArea;
//     }

//     // If we reach this point in the code some part of the face will be swept
//     // during [tSmall, dt-tSmall]. However, it may be the case that there are no
//     // vertex times within the interval. This will happen sometimes for small
//     // time steps where both the initial and the final face-interface
//     // intersection line (FIIL) will be along the same two edges.

//     // Face-interface intersection line (FIIL) to be swept across face
//     DynamicList<point> FIIL(3);
//     // Submerged area at beginning of each sub time interval time
//     scalar initialArea = 0.0;
//     // Running time keeper variable for the integration process
//     scalar time = 0.0;

//     // Special treatment of first sub time interval
//     if (firstTime > 0.0)
//     {
//         // If firstTime > 0 the face is uncut in the time interval
//         // [0, firstTime] and hence fully submerged in fluid A or B.
//         // If Un0 > 0 cell is filling up and it must initially be empty.
//         // If Un0 < 0 cell must initially be fully immersed in fluid A (liquid).
//         time = firstTime;
//         initialArea = magSf * neg(Un0);
//         tIntArea = initialArea * time;
//         cutPoints(faceI, time, pTimes, FIIL);
//     }
//     else
//     {
//         // If firstTime <= 0 then face is initially cut and we must
//         // calculate the initial submerged area and FIIL:
//         time = 0.0;
//         // Note: calcSubFace assumes well-defined 2-point FIIL!!!!
//         calcSubFace(fPts, normal, distance);
//         initialArea = mag(subFaceArea());
//         cutPoints(faceI, time, pTimes, FIIL);
//     }

//     // Making sorted array of all vertex times that are between max(0,firstTime)
//     // and dt and further than tSmall from the previous time.
//     DynamicList<scalar> sortedTimes(pTimes.size());
//     {
//         scalar prevTime = time;
//         const scalar tSmall = max(1e-6*dt, 10*SMALL);
//         for (const label oI : order)
//         {
//             const scalar timeI = pTimes[oI];
//             if (timeI > prevTime + tSmall && timeI <= dt)
//             {
//                 sortedTimes.append(timeI);
//                 prevTime = timeI;
//             }
//         }
//     }

//     // Sweeping all quadrilaterals corresponding to the intervals defined above
//     for (const scalar newTime : sortedTimes)
//     {
//         // New face-interface intersection line
//         DynamicList<point> newFIIL(3);
//         cutPoints(faceI, newTime, pTimes, newFIIL);

//         // quadrilateral area coefficients
//         scalar alpha = 0, beta = 0;

//         quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
//         // Integration of area(t) = A*t^2+B*t from t = 0 to 1
//         tIntArea +=
//             (newTime - time)
//           * (initialArea + sign(Un0)
//           * (alpha / 3.0 + 0.5 * beta));
//         // Adding quad area to submerged area
//         initialArea += sign(Un0) * (alpha + beta);

//         FIIL = newFIIL;
//         time = newTime;
//     }

//     if (lastTime > dt)
//     {
//         // FIIL will end up cutting the face at dt
//         // New face-interface intersection line
//         DynamicList<point> newFIIL(3);
//         cutPoints(faceI, dt, pTimes, newFIIL);

//         // quadrilateral area coefficients
//         scalar alpha = 0, beta = 0;
//         quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
//         // Integration of area(t) = A*t^2+B*t from t = 0 to 1
//         tIntArea +=
//             (dt - time)
//           * (initialArea + sign(Un0) * (alpha / 3.0 + 0.5 * beta));
//     }
//     else
//     {
//         // Interface will leave the face at lastTime and face will be fully
//         // in fluid A or fluid B in the time interval from lastTime to dt.
//         tIntArea += magSf * (dt - lastTime) * pos0(Un0);
//     }

//     return tIntArea;
// }


void Foam::geometricVofExt::SimPLIC::cutFace::quadAreaCoeffs
(
    const DynamicList<point>& pf0,
    const DynamicList<point>& pf1,
    scalar& alpha,
    scalar& beta
) const
{
    // Number of points in provided face-interface intersection lines
    const label np0 = pf0.size();
    const label np1 = pf1.size();

    // quad area coeffs such that area(t) = alpha*t^2 + beta*t.
    // With time interval normalised, we have full quadArea = alpha + beta
    // and time integrated quad area = alpha/3 + beta/2;
    alpha = 0.0;
    beta = 0.0;

    if (np0 && np1)
    {
        // Initialising quadrilateral vertices A, B, C and D
        vector A(pf0[0]);
        vector C(pf1[0]);
        vector B(pf0[0]);
        vector D(pf1[0]);

        if (np0 == 2)
        {
            B = pf0[1];
        }
        else if (np0 > 2)
        {
            WarningInFunction << "Vertex face was cut at pf0 = " << pf0 << endl;
        }

        if (np1 == 2)
        {
            D = pf1[1];
        }
        else if (np1 > 2)
        {
            WarningInFunction << "Vertex face was cut at pf1 = " << pf1 << endl;
        }

        // Swapping pf1 points if pf0 and pf1 point in same general direction
        // (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if (((B - A) & (D - C)) > 0)
        {
            vector tmp = D;
            D = C;
            C = tmp;
        }

        // Defining local coordinates (xhat, yhat) for area integration of swept
        // quadrilateral ABCD such that A = (0,0), B = (Bx,0), C = (Cx,Cy) and
        // D = (Dx,Dy) with Cy = 0 and Dy > 0.

        const scalar Bx = mag(B - A);

        vector xhat(Zero);
        if (Bx > 10 * SMALL)
        {
            // If |AB| > 0 ABCD we use AB to define xhat
            xhat = (B - A) / mag(B - A);
        }
        else if (mag(C - D) > 10 * SMALL)
        {
            // If |AB| ~ 0 ABCD is a triangle ACD and we use CD for xhat
            xhat = (C - D) / mag(C - D);
        }
        else
        {
            return;
        }

        // Defining vertical axis in local coordinates
        vector yhat = D - A;
        yhat -= ((yhat & xhat) * xhat);

        if (mag(yhat) > 10 * SMALL)
        {
            yhat /= mag(yhat);

            const scalar Cx = (C - A) & xhat;
            const scalar Cy = mag((C - A) & yhat);
            const scalar Dx = (D - A) & xhat;
            const scalar Dy = mag((D - A) & yhat);

            // area = ((Cx - Bx)*Dy - Dx*Cy)/6.0 + 0.25*Bx*(Dy + Cy);
            alpha = 0.5 * ((Cx - Bx) * Dy - Dx * Cy);
            beta = 0.5 * Bx * (Dy + Cy);
        }
    }
    else
    {
        WarningInFunction
            << "Vertex face was cut at " << pf0 << " and at "
            << pf1 << endl;
    }
}


void Foam::geometricVofExt::SimPLIC::cutFace::cutPoints
(
    const label faceI,
    const scalar f0,
    const scalarField& pTimes,
    DynamicList<point>& cutPoints
)
{
    const face& f = mesh_.faces()[faceI];
    const label nPoints = f.size();
    scalar f1(pTimes[0]);

    // Snapping vertex value to f0 if very close (needed for 2D cases)
    if (mag(f1 - f0) < 10 * SMALL)
    {
        f1 = f0;
    }

    forAll(f, pi)
    {
        label pi2 = (pi + 1) % nPoints;
        scalar f2 = pTimes[pi2];

        // Snapping vertex value
        if (mag(f2 - f0) < 10 * SMALL)
        {
            f2 = f0;
        }

        if ((f1 < f0 && f2 > f0) || (f1 > f0 && f2 < f0))
        {
            const scalar s = (f0 - f1) / (f2 - f1);
            cutPoints.append
            (
                mesh_.points()[f[pi]]
              + s*(mesh_.points()[f[pi2]] - mesh_.points()[f[pi]])
            );
        }
        else if (f1 == f0)
        {
            cutPoints.append(mesh_.points()[f[pi]]);
        }
        f1 = f2;
    }

    if (cutPoints.size() > 2)
    {
        WarningInFunction
            << "cutPoints = " << cutPoints
            << " for pts = " << f.points(mesh_.points())
            << ", f - f0 = " << f - f0 << " and f0 = " << f0
            << endl;
    }
}


void Foam::geometricVofExt::SimPLIC::cutFace::cutPoints
(
    const pointField& pts,
    const scalarField& f,
    const scalar f0,
    DynamicList<point>& cutPoints
)
{
    const label nPoints = pts.size();
    scalar f1(f[0]);

    // Snapping vertex value to f0 if very close (needed for 2D cases)
    if (mag(f1 - f0) < 10 * SMALL)
    {
        f1 = f0;
    }

    forAll(pts, pi)
    {
        label pi2 = (pi + 1) % nPoints;
        scalar f2 = f[pi2];

        // Snapping vertex value
        if (mag(f2 - f0) < 10 * SMALL)
        {
            f2 = f0;
        }

        if ((f1 < f0 && f2 > f0) || (f1 > f0 && f2 < f0))
        {
            const scalar s = (f0 - f1) / (f2 - f1);
            cutPoints.append(pts[pi] + s * (pts[pi2] - pts[pi]));
        }
        else if (f1 == f0)
        {
            cutPoints.append(pts[pi]);
        }
        f1 = f2;
    }

    if (cutPoints.size() > 2)
    {
        WarningInFunction
            << "cutPoints = " << cutPoints << " for pts = " << pts
            << ", f - f0 = " << f - f0 << " and f0 = " << f0
            << endl;
    }
}

// ************************************************************************* //