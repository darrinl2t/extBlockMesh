/*---------------------------------------------------------------------------*\
  extBlockMesh
  Copyright (C) 2014 Etudes-NG
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

#include "SmootherControl.H"

#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SmootherControl::SmootherControl
(
    const dictionary& smootherDict,
    const word& subDictName
)
{
    // Get smoother parameters
    const dictionary& ctrls = smootherDict.subDict(subDictName);

    ctrls.readEntry("maxIterations", _maxIterations);
    ctrls.readEntry("transformParameter", _transformParam);
    ctrls.readEntry("meanImprovTol", _meanImprovTol);
    ctrls.readEntry("maxMinCycleNoChange", _maxMinCycleNoChange);
    ctrls.readEntry("meanRelaxationTable", _meanRelaxTable);
    ctrls.readEntry("minRelaxationTable", _minRelaxTable);
    ctrls.readEntry("snapRelaxationTable", _snapRelaxTable);
    ctrls.readEntry("ratioWorstQualityForMin", _ratioForMin);

    if (_meanRelaxTable.last() > VSMALL)
    {
        _meanRelaxTable.append(0);
    }
    if (_minRelaxTable.last() > VSMALL)
    {
        _minRelaxTable.append(0);
    }

    if (_snapRelaxTable.last() > VSMALL)
    {
        _snapRelaxTable.append(0);
    }
//    if (_snapRelaxTable.first() < (1.0 - VSMALL))
//    {
//        scalarList relax(_snapRelaxTable.size() + 1);
//        relax[0] = 1.0;

//        forAll(_snapRelaxTable, relaxI)
//        {
//            relax[relaxI + 1] = _snapRelaxTable[relaxI];
//        }
//        _snapRelaxTable = relax;
//    }

    Info<< nl
        << "  smoothControls:"  << nl
        << "    - Max iterations             : " << _maxIterations  << nl
        << "    - Tranformation parameter    : " << _transformParam << nl
        << "    - Mean improvement tolerance : " << _meanImprovTol << nl
        << "    - Max ineffective iteration  : " << _maxMinCycleNoChange << nl
        << "    - Mean relaxation table      : " << _meanRelaxTable << nl
        << "    - Min relaxation table       : " << _minRelaxTable << nl
        << "    - Snap relaxation table      : " << _snapRelaxTable << nl
        << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
