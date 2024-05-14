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

#include "cutCell.H"
#include "mergePoints.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::geometricVofExt::SimPLIC::cutCell::calcInterfaceCentreAndArea()
{
    point fCentre(point::zero);
    label nEdgePoints(0);

    for (const DynamicList<point>& edgePoints : interfaceEdges_)
    {
        for (const point& p : edgePoints)
        {
            fCentre += p;
            nEdgePoints++;
        }
    }

    if (nEdgePoints > 0)
    {
        fCentre /= nEdgePoints;
    }

    vector sumN(vector::zero);
    scalar sumA(0.0);
    vector sumAc(vector::zero);

    forAll(interfaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints(interfaceEdges_[ei]);
        const label nPoints(edgePoints.size());
        for (label pi = 0; pi < nPoints - 1; pi++)
        {
            const point& nextPoint(edgePoints[pi + 1]);

            vector c(edgePoints[pi] + nextPoint + fCentre);
            vector n
            (
                (nextPoint - edgePoints[pi]) ^ (fCentre - edgePoints[pi])
            );
            scalar a(mag(n));

            // Edges may have different orientation
            sumN += Foam::sign(n & sumN) *  n;
            sumA += a;
            sumAc += a * c;
        }
    }

    // This is to deal with zero-area faces. Mark very small faces
    // to be detected in e.g., processorPolyPatch.
    if (sumA < ROOTVSMALL)
    {
        interfaceCentre_ = fCentre;
        interfaceArea_ = vector::zero;
    }
    else
    {
        interfaceCentre_ = (1.0 / 3.0) * sumAc / sumA;
        interfaceArea_ = 0.5 * sumN;
    }

    // Check faceArea direction and change if not pointing in gas/secondary
    if ((interfaceArea_ & (interfaceCentre_ - subCellCentre_)) < 0.0)
    {
        interfaceArea_ *= (-1.0);
    }
}


void Foam::geometricVofExt::SimPLIC::cutCell::calcSubCellCentreAndVolume()
{
    // Clear the fields for accumulation
    subCellCentre_ = point::zero;
    subCellVolume_ = 0.0;

    // Estimate the approximate cell centre as the average of face centres
    vector cEst(average(cutFaceCentres_));

    // Contribution to subcell centre and volume from cut faces
    forAll(cutFaceCentres_, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol
        (
            max
            (
                mag(cutFaceAreas_[facei] & (cutFaceCentres_[facei] - cEst)),
                VSMALL
            )
        );

        // Calculate face-pyramid centre
        vector pc(0.75 * cutFaceCentres_[facei] + 0.25 * cEst);

        // Accumulate volume-weighted face-pyramid centre
        subCellCentre_ += pyr3Vol * pc;

        // Accumulate face-pyramid volume
        subCellVolume_ += pyr3Vol;
    }

    subCellCentre_ /= subCellVolume_;
    subCellVolume_ /= 3.0; // formula of pyramid
}


void Foam::geometricVofExt::SimPLIC::cutCell::getLocalPointFieldAndFaceList
(
    bool splitWarpedFace
)
{
    localPoints_.clear();
    localFaces_.clear();

    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    // Refs to mesh points, faces, cells
    const pointField& points(mesh_.points());
    const faceList& faces(mesh_.faces());
    const cellList& cells(mesh_.cells());
    const surfaceVectorField& Cf(mesh_.Cf());

    const vectorField& fCtrs(mesh_.faceCentres());
    const labelList& own(mesh_.faceOwner());

    // Cell cellI_ info
    const cell& c(cells[cellI_]);
    const labelList globalPointLabels(c.labels(faces));

    localPoints_.clear();
    localFaces_.clear();

    forAll(globalPointLabels, pointi)
    {
        localPoints_.append(points[globalPointLabels[pointi]]);
    }

    if (splitWarpedFace)
    {
        forAll(c, facei)
        {
            if (faceFlatness_[c[facei]] > (1.0-TSMALL)) // facei is flat
            {
                const face& fa(faces[c[facei]]);
                face localFa(fa.size());
                localFa.clear();

                forAll(fa, pointi)
                {
                    localFa.append(globalPointLabels.find(fa[pointi]));
                }

                localFaces_.append
                (
                    cellI_ == own[c[facei]] ? localFa : localFa.reverseFace()
                );
            }
            else    // face is not flat, triangular decomposition
            {
                const face& fa(faces[c[facei]]);

                localPoints_.append(fCtrs[c[facei]]);

                forAll(fa, pointi)
                {
                    const label nextPointi((pointi + 1) % fa.size());

                    face localFa(3);
                    //localFa.clear();

                    localFa[0] = localPoints_.size()-1;
                    localFa[1] = globalPointLabels.find(fa[pointi]);
                    localFa[2] = globalPointLabels.find(fa[nextPointi]);

                    localFaces_.append
                    (
                        cellI_ == own[c[facei]] ? localFa : localFa.reverseFace()
                    );
                }
            }
        }
    }
    else
    {
        forAll(c, facei)
        {
            const face& fa(faces[c[facei]]);
            face localFa(fa.size());
            localFa.clear();

            forAll(fa, pointi)
            {
                localFa.append(globalPointLabels.find(fa[pointi]));
            }

            localFaces_.append
            (
                cellI_ == own[c[facei]] ? localFa : localFa.reverseFace()
            );
        }
    }
}


void Foam::geometricVofExt::SimPLIC::cutCell::updateSubCellPointsandFaces()
{
    // Append interface points
    const DynamicList<point>& interfacePts(interfacePoints());
    if (interfacePoints_.size() > 0)
    {
        subCellFaces_.append
        (
            face
            (
                identity
                (
                    interfacePoints_.size(),
                    subCellPoints_.size()
                )
            )
        );

        subCellPoints_.append(interfacePoints_);
    }

    // Merge duplicated points and update point labels of faces
    labelList pointToUnique;
    inplaceMergePoints<DynamicList<point>>
    (
        subCellPoints_,
        10.0*SMALL,
        false,
        pointToUnique
    );

    forAll (subCellFaces_, i)
    {
        face& facei(subCellFaces_[i]);

        forAll (facei, pti)
        {
            facei[pti] = pointToUnique[facei[pti]];
        }

        // Correct face orientations
        const point faceiC(facei.centre(subCellPoints_));
        const vector faceiN(facei.normal(subCellPoints_));
        const vector cellCToFaceiC(faceiC - subCellCentre_);
        if ((cellCToFaceiC & faceiN) < 0.0)
        {
            facei = facei.reverseFace();
        }
    }

    isSubCellPointsAndFacesUpdated_ = true;
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::geometricVofExt::SimPLIC::cutCell::cutCell
(
    const fvMesh& mesh,
    const scalarField& faceFlatness,
    const volScalarField& alpha1,
    const volVectorField& interfaceN,
    volScalarField& interfaceD,
    volVectorField& interfaceC,
    volVectorField& interfaceS
)
:
    mesh_(mesh),
    faceFlatness_(faceFlatness),
    alpha1_(alpha1),
    interfaceN_(interfaceN),
    interfaceD_(interfaceD),
    interfaceC_(interfaceC),
    interfaceS_(interfaceS),
    cutFace_(mesh, faceFlatness),
    cellStatus_(-1),
    cutFaceCentres_(10),
    cutFaceAreas_(10),
    subCellPoints_(50),
    subCellFaces_(30),
    isSubCellPointsAndFacesUpdated_(false),
    subCellCentre_(point::zero),
    subCellVolume_(0.0),
    VOF_(0.0),
    cellI_(-1),
    interfaceEdges_(10),
    interfacePoints_(10),
    interfaceCentre_(point::zero),
    interfaceArea_(vector::zero),
    localPoints_(100),
    localFaces_(100)
{
    mesh.C();
    mesh.V();
    mesh.Cf();
    mesh.magSf();

    clearStorage();
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //


Foam::label Foam::geometricVofExt::SimPLIC::cutCell::calcSubCell
(
    const label cellI,
    const vector& normal,
    const scalar distance,
    bool splitWarpedFace
)
{
    clearStorage();

    cellI_ = cellI;

    const scalar TSMALL(10.0*SMALL);    // Tolerance

    const vectorField& C(mesh_.cellCentres());
    const scalarField& V(mesh_.cellVolumes());
    const pointField& points(mesh_.points());

    bool fullySubmerged(true);
    bool fullyEmpty(true);

    label nSubmergedFaces(0);

    if (splitWarpedFace)
    {
        forAll(localFaces_, i)
        {
            const pointField fPts(localFaces_[i].points(localPoints_));
            const label faceStatus
            (
                cutFace_.calcSubFace(fPts, normal, distance)
            );

            if (faceStatus == 0)    // face is cut
            {
                cutFaceCentres_.append(cutFace_.subFaceCentre());
                cutFaceAreas_.append(cutFace_.subFaceArea());
                interfaceEdges_.append(cutFace_.interfacePoints());

                subCellFaces_.append
                (
                    face
                    (
                        identity
                        (
                            cutFace_.subFacePoints().size(),
                            subCellPoints_.size()
                        )
                    )
                );
                subCellPoints_.append(cutFace_.subFacePoints());

                fullySubmerged = false;
                fullyEmpty = false;
            }
            else if (faceStatus == -1) // face is fully submerged
            {
                cutFaceCentres_.append(cutFace_.subFaceCentre());
                cutFaceAreas_.append(cutFace_.subFaceArea());

                subCellFaces_.append
                (
                    face
                    (
                        identity
                        (
                            fPts.size(),
                            subCellPoints_.size()
                        )
                    )
                );
                subCellPoints_.append(fPts);

                fullyEmpty = false;
                nSubmergedFaces++;
            }
            else
            {
                fullySubmerged = false;
            }
        }
    }
    else
    {
        const cell& c(mesh_.cells()[cellI]);

        forAll(c, fi)
        {
            const label faceI(c[fi]);
            const label faceStatus
            (
                cutFace_.calcSubFace(faceI, normal, distance)
            );

            if (faceStatus == 0)    // face is cut
            {
                cutFaceCentres_.append(cutFace_.subFaceCentre());
                cutFaceAreas_.append(cutFace_.subFaceArea());
                interfaceEdges_.append(cutFace_.interfacePoints());

                subCellFaces_.append
                (
                    face
                    (
                        identity
                        (
                            cutFace_.subFacePoints().size(),
                            subCellPoints_.size()
                        )
                    )
                );
                subCellPoints_.append(cutFace_.subFacePoints());

                fullySubmerged = false;
                fullyEmpty = false;
            }
            else if (faceStatus == -1) // face is fully submerged
            {
                cutFaceCentres_.append(cutFace_.subFaceCentre());
                cutFaceAreas_.append(cutFace_.subFaceArea());

                subCellFaces_.append
                (
                    face
                    (
                        identity
                        (
                            mesh_.faces()[faceI].size(),
                            subCellPoints_.size()
                        )
                    )
                );
                subCellPoints_.append(mesh_.faces()[faceI].points(points));

                fullyEmpty = false;
                nSubmergedFaces++;
            }
            else
            {
                fullySubmerged = false;
            }
        }
    }

    if (!fullySubmerged && !fullyEmpty) // cell cut at least at one face
    {
        cellStatus_ = 0;

        calcInterfaceCentreAndArea();

        // In the rare but occuring cases where a cell is only touched at a
        // point or a line the interfaceArea_ will have zero length and
        // here the cell should be treated as either completely empty or full.
        if (mag(interfaceArea_) < TSMALL)
        {
            if (nSubmergedFaces == 0)   // Fully empty cell
            {
                cellStatus_ = 1;
                subCellCentre_ = point::zero;
                subCellVolume_ = 0.0;
                VOF_ = 0.0;

                return cellStatus_;
            }
            else    // Fully submerged cell
            {
                cellStatus_ = -1;
                subCellCentre_ = C[cellI];
                subCellVolume_ = V[cellI];
                VOF_ = 1.0;

                return cellStatus_;
            }
        }

        cutFaceCentres_.append(interfaceCentre_);
        cutFaceAreas_.append(interfaceArea_);

        // Calc volume and sub cell centre
        calcSubCellCentreAndVolume();

        VOF_ = subCellVolume_ / V[cellI];
    }
    else if (fullyEmpty)  // Fully empty cell
    {
        cellStatus_ = 1;
        subCellCentre_ = point::zero;
        subCellVolume_ = 0.0;
        VOF_ = 0.0;
    }
    else if (fullySubmerged)    // Fully submerged cell
    {
        cellStatus_ = -1;
        subCellCentre_ = C[cellI];
        subCellVolume_ = V[cellI];
        VOF_ = 1.0;
    }

    return cellStatus_;
}


const Foam::DynamicList<Foam::point>&
Foam::geometricVofExt::SimPLIC::cutCell::interfacePoints()
{
    if (interfacePoints_.empty())
    {
        if (cellStatus_ == 0)
        {
            const label nEdges(interfaceEdges_.size());

            // Define local coordinates with zhat along interface normal
            // and xhat from interface centre to first point in interfaceEdges_
            const vector zhat(interfaceArea_ / mag(interfaceArea_));
            vector xhat(interfaceEdges_[0][0] - interfaceCentre_);
            xhat = (xhat - (xhat & zhat) * zhat);
            xhat.normalise();
            vector yhat(zhat ^ xhat);
            yhat.normalise();

            // Calculate interface point angles in local coordinates
            DynamicList<point> unsortedPlicFacePoints(3 * nEdges);
            DynamicList<scalar> unsortedPlicFacePointAngles(3 * nEdges);
            for (const DynamicList<point>& edgePoints : interfaceEdges_)
            {
                for (const point& p : edgePoints)
                {
                    unsortedPlicFacePoints.append(p);
                    unsortedPlicFacePointAngles.append
                    (
                        Foam::atan2
                        (
                            ((p - interfaceCentre_) & yhat),
                            ((p - interfaceCentre_) & xhat)
                        )
                    );
                }
            }

            // Sort interface points by angle and insert into interfacePoints_
            labelList order(unsortedPlicFacePointAngles.size());
            Foam::sortedOrder(unsortedPlicFacePointAngles, order);
            interfacePoints_.append(unsortedPlicFacePoints[order[0]]);
            for (label pi = 1; pi < order.size(); pi++)
            {
                if
                (
                    mag
                    (
                        unsortedPlicFacePointAngles[order[pi]]
                      - unsortedPlicFacePointAngles[order[pi-1]]
                    ) > 1e-8
                )
                {
                    interfacePoints_.append(unsortedPlicFacePoints[order[pi]]);
                }
            }
        }
        else
        {
            interfacePoints_.clear();
        }
    }

    return interfacePoints_;
}


Foam::label Foam::geometricVofExt::SimPLIC::cutCell::findSignedDistance
(
    const label cellI,
    bool splitWarpedFace
)
{
    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    // Cell cellI_ info
    cellI_ = cellI;
    const scalar alphaI(alpha1_.primitiveField()[cellI_]);
    const vector& normalI(interfaceN_.primitiveField()[cellI_]);

    if (mag(normalI) < TSMALL)
    {
        return sign(0.5 - alphaI);
    }

    const vectorField& C(mesh_.cellCentres());
    const scalarField& V(mesh_.cellVolumes());
    const pointField& points(mesh_.points());
    // const faceList& faces(mesh_.faces());
    // const cellList& cells(mesh_.cells());
    // const cell& c(cells[cellI_]);

    label nPoints(mesh_.cellPoints(cellI).size());
    if (splitWarpedFace)
    {
        getLocalPointFieldAndFaceList(true);
        nPoints = localPoints_.size();
    }

    // Calculate and sort vertex distances
    scalarField vertexDistances(nPoints);

    if (splitWarpedFace)
    {
        forAll(localPoints_, pointi)
        {
            vertexDistances[pointi] = -(normalI & localPoints_[pointi]);
        }
    }
    else
    {
        const labelList& pLabels(mesh_.cellPoints(cellI));

        forAll(pLabels, pointi)
        {
            vertexDistances[pointi] = -(normalI & points[pLabels[pointi]]);
        }
    }

    labelList vertexDistanceOrder(vertexDistances.size());
    Foam::sortedOrder
    (
        vertexDistances,
        vertexDistanceOrder,
        typename UList<scalar>::greater(vertexDistances)
    );

    // Binary Bracketing (BB) procedure
    //   Details can be found in:
    //       \verbatim
    //           L{\'o}pez, Joaqu{\'\i}n and Hern{\'a}ndez, J. (2008).
    //           Analytical and geometrical tools for 3D volume of fluid
    //           methods in general grids
    //           Journal of Computational Physics
    //           doi 10.1016/j.jcp.2008.03.010
    //           url https://doi.org/10.1016/j.jcp.2008.03.010
    //       \endverbatim
    scalar lowDistance(vertexDistances[vertexDistanceOrder.first()]);
    scalar upDistance(vertexDistances[vertexDistanceOrder.last()]);
    label lowLabel(0);
    label upLabel(vertexDistances.size() - 1);
    scalar lowAlpha(0.0);
    scalar upAlpha(1.0);

    scalar midDistance, midLabel, midAlpha;

    while ((upLabel - lowLabel) > 1)
    {
        midLabel = round(0.5 * (upLabel + lowLabel));
        midDistance = vertexDistances[vertexDistanceOrder[midLabel]];
        calcSubCell(cellI, normalI, midDistance, splitWarpedFace);
        midAlpha = volumeOfFluid();

        if (mag(midAlpha - alphaI) < TSMALL)
        {
            interfaceD_.primitiveFieldRef()[cellI_] = midDistance;
            interfaceC_.primitiveFieldRef()[cellI_] = interfaceCentre_;
            interfaceS_.primitiveFieldRef()[cellI_] = interfaceArea_;

            return cellStatus_;
        }

        if (midAlpha > alphaI)
        {
            upLabel = midLabel;
            upDistance = midDistance;
            upAlpha = midAlpha;
        }
        else
        {
            lowLabel = midLabel;
            lowDistance = midDistance;
            lowAlpha = midAlpha;
        }
    }

    if (mag(lowDistance - upDistance) < TSMALL)
    {
        const scalar midDistance(0.5 * (lowDistance + upDistance));
        calcSubCell(cellI, normalI, midDistance, splitWarpedFace);
        interfaceD_.primitiveFieldRef()[cellI_] = midDistance;
        interfaceC_.primitiveFieldRef()[cellI_] = interfaceCentre_;
        interfaceS_.primitiveFieldRef()[cellI_] = interfaceArea_;

        return cellStatus_;
    }

    // Standard distance-finding algorithm
    //   Details can be found in:
    //       \verbatim
    //           Dai, Dezhi and Tong, Albert Y. (2019).
    //           Analytical interface reconstruction algorithms in the PLIC‐VOF
    //           method for 3D polyhedral unstructured meshes
    //           International Journal for Numerical Methods in Fluids
    //           doi 10.1002/fld.4750
    //           url https://doi.org/10.1002/fld.4750
    //       \endverbatim

    // Fraction value of the prismatoid
    const scalar alphaPrismatoid(upAlpha - lowAlpha);

    // Finding 2 additional points
    const scalar deltaDistance((upDistance - lowDistance) / 3.0);

    const scalar distanceOneThird(lowDistance + deltaDistance);
    calcSubCell(cellI, normalI, distanceOneThird, splitWarpedFace);
    const scalar alphaOneThird(volumeOfFluid() - lowAlpha);



    const scalar distanceTwoThirds(lowDistance + 2.0 * deltaDistance);
    calcSubCell(cellI, normalI, distanceTwoThirds, splitWarpedFace);
    const scalar alphaTwoThirds(volumeOfFluid() - lowAlpha);



    // Calculate a, b c and d by using Eqs. (20) & (17)
    scalar a, b, c, d;

    a =  13.5 * alphaOneThird - 13.5 * alphaTwoThirds + 4.5 * alphaPrismatoid;

    b = -22.5 * alphaOneThird + 18.0 * alphaTwoThirds - 4.5 * alphaPrismatoid;

    c =   9.0 * alphaOneThird -  4.5 * alphaTwoThirds + 1.0 * alphaPrismatoid;

    d = lowAlpha - alphaI;

    // Find the root of Eq. (16) by using Newton method, max iteration is 100
    scalar lambda(0.5);    // Initial guess is set to 0.5
    for (label iter = 0; iter < 100; iter++)
    {
        const scalar func(a * pow3(lambda) + b * sqr(lambda) + c * lambda + d);
        const scalar funcPrime(3.0 * a * sqr(lambda) + 2.0 * b * lambda + c);
        const scalar lambdaNew(lambda - (func / funcPrime));

        if (mag(lambdaNew - lambda) < TSMALL) // Convergence tolerance is 1e-14
        {
            break;
        }

        lambda = lambdaNew;
    }

    // Calculate signed distance by using Eq. (18)
    const scalar distance0(lowDistance - lambda * (lowDistance - upDistance));

    // Update
    calcSubCell(cellI, normalI, distance0, splitWarpedFace);

    interfaceD_.primitiveFieldRef()[cellI_] = distance0;
    interfaceC_.primitiveFieldRef()[cellI_] = interfaceCentre_;
    interfaceS_.primitiveFieldRef()[cellI_] = interfaceArea_;

    return cellStatus_;
}


// ************************************************************************* //