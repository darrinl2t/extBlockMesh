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

#include "SmootherBoundaryLayer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherBoundaryLayer::SmootherBoundaryLayer()
:
    _nbLayers(0),
    _expansionRatio(1),
    _finalLayerThickness(1),
    _relativeSize(true)
{}


Foam::SmootherBoundaryLayer::SmootherBoundaryLayer(const dictionary& dict)
:
    SmootherBoundaryLayer()
{
    dict.readEntry("nSurfaceLayers", _nbLayers);
    dict.readEntry("expansionRatio", _expansionRatio);
    dict.readEntry("finalLayerThickness", _finalLayerThickness);
    dict.readEntry("relativeSizes", _relativeSize);

    Info<< "        - Number of BL       : " << _nbLayers << nl
        << "        - Expansion ratio    : " << _expansionRatio << nl
        << "        - Relative size      : " << _relativeSize << nl
        << "        - Final thickness    : " << _finalLayerThickness << nl;
}


// ************************************************************************* //
