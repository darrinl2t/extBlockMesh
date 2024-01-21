/*---------------------------------------------------------------------------*\
  extBlockMesh
  Copyright (C) 2014 Etudes-NG
  Copyright (C) 2020 OpenCFD Ltd.
  ---------------------------------
License
    This file is part of extBlockMesh.

    extBlockMesh is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    extBlockMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with extBlockMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "MeshSmoother.H"
#include "SmootherPoint.H"
#include "SmootherCell.H"
#include "SmootherControl.H"
#include "SmootherParameter.H"
#include "SmootherBoundary.H"

#include "Time.H"
#include "IOmanip.H"

#include "polyMeshFields.H"

#include <algorithm>
#include <cmath>
#include <memory>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Write mesh with error check
    static inline void writeMesh(polyMesh* mesh)
    {
        if (mesh && !mesh->write())
        {
            FatalErrorInFunction
                << "Failed writing polyMesh."
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //

void Foam::MeshSmoother::analyseMeshQuality()
{
    label celli = 0;
    for (const cellShape& shape : _polyMesh->cellShapes())
    {
        _cellQuality[celli] = SmootherCell::quality(shape);
        ++celli;
    }
}


void Foam::MeshSmoother::analyseMeshQuality(const labelHashSet &cell)
{
    forAllConstIter(labelHashSet, cell, cellIter)
    {
        const label celli = cellIter.key();
        const cellShape& shape = _polyMesh->cellShapes()[celli];

        _cellQuality[celli] = SmootherCell::quality(shape);
    }
}


void Foam::MeshSmoother::qualityStats()
{
    scalar minQuality = 1.0;
    scalar meanQuality = 0.0;

    scalarList pQSum(_polyMesh->nPoints(), Zero);
    forAll(_polyMesh->cells(), cellI)
    {
        const scalar cQ = _cellQuality[cellI];
        if (cQ < minQuality)
        {
            minQuality = cQ;
        }

        meanQuality += cQ;

        forAll(_polyMesh->cellPoints()[cellI], pointI)
        {
            pQSum[_polyMesh->cellPoints()[cellI][pointI]] += cQ;
        }
    }

    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->setQuality
        (
            pQSum[ptI]/_polyMesh->pointCells()[ptI].size()
        );
    }

    meanQuality /= _polyMesh->nCells();

    _param->setMinQual(minQuality);
    _param->setMeanQual(meanQuality);
}


Foam::labelHashSet Foam::MeshSmoother::addTransformedElementNodeWeight()
{
    labelHashSet transformedPoints;
    forAll(_polyMesh->cells(), cellI)
    {
        const scalar cQ = _cellQuality[cellI];

        if (cQ <= _param->transformationTreshold())
        {
            if (cQ < VSMALL)
            {
                FatalErrorInFunction
                    << "Quality of cell " << cellI << " is null" << nl
                    << exit(FatalError);
            }

            const cellShape& cS = _polyMesh->cellShapes()[cellI];

            const pointField newCellPoints =
                SmootherCell::geometricTransform(cS);

            forAll(cS, pointI)
            {
                SmootherPoint* pt = _bnd->pt(cS[pointI]);

                // compute the associated weight
                const label nNei = _polyMesh->pointCells(cS[pointI]).size();
                const scalar weight = std::sqrt(pt->avgQual()/(nNei*cQ));

                // add the weight to temporary weighted sum
                pt->addWeight(weight, newCellPoints[pointI]);

                // add point to set of tranformed points
                transformedPoints.insert(cS[pointI]);
            }
        }
    }
    return transformedPoints;
}

void Foam::MeshSmoother::addUnTransformedElementNodeWeight(labelHashSet &tp)
{
    // Add Untransformed Element Nodes And Weights
    forAll (_polyMesh->cells(), cellI)
    {
        if (untransformedAndhavePointTransformed(cellI, tp))
        {
            const scalar cQ = _cellQuality[cellI];

            if (cQ < VSMALL)
            {
                FatalErrorInFunction
                    << "Quality of cell " << cellI << " is null" << nl
                    << exit(FatalError);
            }

            const cellShape& cS = _polyMesh->cellShapes()[cellI];

            forAll(cS, pointI)
            {
                SmootherPoint* pt = _bnd->pt(cS[pointI]);

                // compute the associated weight
                const label nNei = _polyMesh->pointCells(cS[pointI]).size();
                const scalar weight = std::sqrt(pt->avgQual()/(nNei*cQ));

                // add the weight to temporary weighted sum
                pt->addWeight(weight);
            }
        }
    }
}

bool Foam::MeshSmoother::untransformedAndhavePointTransformed
(
    const label cellI,
    const labelHashSet& tp
)
{
    if (_cellQuality[cellI] > _param->transformationTreshold())
    {
        const cellShape& cS = _polyMesh->cellShapes()[cellI];
        forAll(cS, ptI)
        {
            if (tp.found(cS[ptI]))
            {
                // Have point transformed
                return true;
            }
        }
    }

    return false;
}

void Foam::MeshSmoother::iterativeNodeRelaxation
(
    labelHashSet &tP,
    const scalarList &r
)
{
    // Reset relaxation level
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->resetRelaxationLevel();
    }
    _param->setNbMovedPoints(tP.size());

    label nbRelax = 0;
    while (!tP.empty())
    {
        ++nbRelax;
        // Relax all the points marked for move
        labelHashSet modifiedCells;
        forAllConstIter(labelHashSet, tP, ptI)
        {
            _bnd->pt(ptI.key())->relaxPoint(r);

            forAll(_polyMesh->pointCells()[ptI.key()], cellI)
            {
                modifiedCells.insert(_polyMesh->pointCells()[ptI.key()][cellI]);
            }
        }

        // compute quality with relaxed points
        tP.clearStorage();
        analyseMeshQuality(modifiedCells);
        forAllConstIter(labelHashSet, modifiedCells, cellIter)
        {
            const label celli = cellIter.key();

            if (_cellQuality[celli] < VSMALL)
            {
                tP.insert(_polyMesh->cellPoints()[celli]);
            }
        }

        // Increase the relaxation level for invalid points
        forAllIter(labelHashSet, tP, pointIter)
        {
            const label pointi = pointIter.key();

            _bnd->pt(pointi)->addRelaxLevel(r);
        }
    }
    _param->setNbRelaxations(nbRelax);
}

bool Foam::MeshSmoother::runIteration()
{
    _param->resetUpdateTime();

    // Store mean and min quality before iteration
    const scalar minQ = _param->minQual();
    const scalar meanQ = _param->meanQual();

    if (!_bnd->unSnapedPoints().empty())
    {
        snapSmoothing();
    }
    else
    {
        GETMeSmoothing();
    }

    // Compute new min and avg quality
    qualityStats();
    _param->printStatus(_bnd->unSnapedPoints().size());

    const bool asUnSnaped = _bnd->unSnapedPoints().empty();
    return _param->setSmoothCycle(meanQ, minQ, asUnSnaped, this);
}


Foam::pointField Foam::MeshSmoother::getMovedPoints() const
{
    pointField pt(_polyMesh->nPoints());

    forAll(pt, ptI)
    {
        pt[ptI] = _bnd->pt(ptI)->getRelaxedPoint();
    }
    return pt;
}


void Foam::MeshSmoother::GETMeSmoothing()
{
    //-------------------------------------------------------------------------

    // Reset all points
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->laplaceReset();
    }

    // LaplaceSmooth boundary points
    labelHashSet snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->featLaplaceSmooth();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());

    // Snap boundary points
    snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->snap();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());

    //-------------------------------------------------------------------------

    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->GETMeReset();
    }

    labelHashSet transformedPoints = addTransformedElementNodeWeight();

    if (transformedPoints.size() != _polyMesh->nPoints())
    {
        addUnTransformedElementNodeWeight(transformedPoints);
    }

    // Compute new point
    forAllConstIter(labelHashSet, transformedPoints, ptI)
    {
        _bnd->pt(ptI.key())->GETMeSmooth();
    }

    iterativeNodeRelaxation(transformedPoints, _param->relaxationTable());

//    _bnd->writeAllSurfaces(_param->getIterNb());
}


void Foam::MeshSmoother::snapSmoothing()
{
    // Reset all points
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->laplaceReset();
    }

    // LaplaceSmooth interior points
    labelHashSet laplacePoints = _bnd->interiorPoints();
    forAllConstIter(labelHashSet, laplacePoints, ptI)
    {
        _bnd->pt(ptI.key())->laplaceSmooth();
    }
    iterativeNodeRelaxation(laplacePoints, _ctrl->snapRelaxTable());

    // Snap boundary points
    labelHashSet snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->snap();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());

    // Remove points from unsnaped point list if snaped
    forAll(_polyMesh->points(), ptI)
    {
        _bnd->pt(ptI)->needSnap();
    }

    // LaplaceSmooth boundary points
    snapPoints = _bnd->featuresPoints();
    forAllConstIter(labelHashSet, snapPoints, ptI)
    {
        _bnd->pt(ptI.key())->featLaplaceSmooth();
    }
    iterativeNodeRelaxation(snapPoints, _ctrl->snapRelaxTable());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MeshSmoother::MeshSmoother
(
    polyMesh *mesh,
    const dictionary& smootherDict
)
:
    _polyMesh(mesh),
    _cellQuality(mesh->nCells(), Zero)
{
    const scalar time = _polyMesh->time().elapsedCpuTime();

    _ctrl = new SmootherControl(smootherDict);
    _param = new SmootherParameter(_ctrl, _polyMesh);

    const dictionary& snapDict = smootherDict.subDict("snapControls");
    _bnd = new SmootherBoundary(snapDict, _polyMesh);

    SmootherCell::setStaticItems(_bnd, _ctrl->transformationParameter());
    SmootherPoint::setStaticItems(_bnd, _param, _polyMesh);

    // Analyse initial quality
    analyseMeshQuality();
    qualityStats();

    //snapFeatures();

    Info<< "Smoother initialized in "
        << _polyMesh->time().elapsedCpuTime() - time << " s" << nl << nl
        << "Smooth the mesh" << nl << nl;

    _param->resetUpdateTime();
    _param->printHeaders();
    _param->setSmoothCycle(!_bnd->unSnapedPoints().empty());
    _param->printStatus(_bnd->unSnapedPoints().size());
    _param->setIterNb();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::MeshSmoother::~MeshSmoother()
{
    delete _param;
    delete _bnd;
    delete _ctrl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MeshSmoother::update()
{
    while (runIteration()) {}

    // Update the mesh with new points
    _polyMesh->movePoints(getMovedPoints());

    _param->printStats();
}


Foam::label Foam::MeshSmoother::updateAndWrite
(
    Time& runTime,
    const bool withQuality
)
{
    std::unique_ptr<polyMeshScalarField> qualityFieldPtr;

    if (withQuality)
    {
        qualityFieldPtr.reset
        (
            new polyMeshScalarField
            (
                IOobject
                (
                    "meshQuality",
                    runTime.timeName(),
                    *_polyMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                *_polyMesh,
                dimless,
                _cellQuality  // Copy in current quality field
            )
        );
    }


    label nWritten = 1;
    writeMesh(_polyMesh);

    while (runIteration())
    {
        // Update the mesh with new points
        _polyMesh->movePoints(getMovedPoints());

        ++runTime;
        ++nWritten;

        if (qualityFieldPtr)
        {
            // Update quality for output
            qualityFieldPtr->field() = _cellQuality;
        }

        writeMesh(_polyMesh);
    }

    _param->printStats();

    return nWritten;
}


Foam::scalar Foam::MeshSmoother::getTransformationTreshold() const
{
    labelList order;
    Foam::sortedOrder(_cellQuality, order);

    const label minIndex = std::floor(order.size() * _ctrl->ratioForMin());

    return _cellQuality[order[minIndex]];
}


// ************************************************************************* //
