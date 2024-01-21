/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2014 Etudes-NG
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    blockMesh

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
    \a constant/polyMesh/blockMeshDict
    (or \a constant/\<region\>/polyMesh/blockMeshDict).

Usage

    - blockMesh [OPTION]

    \param -blockTopology \n
    Write the topology as a set of edges in OBJ format.

    \param -region \<name\> \n
    Specify an alternative mesh region.

    \param -dict \<filename\> \n
    Specify alternative dictionary for the block mesh description.

\*---------------------------------------------------------------------------*/
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

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "attachPolyTopoChanger.H"
#include "emptyPolyPatch.H"

#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"

#include "MeshSmoother.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addBoolOption
    (
        "writeSteps",
        "Write mesh at different smoothing steps"
    );
    argList::addBoolOption
    (
        "write-quality",
        "Write mesh and quality field at different smoothing steps"
    );
    argList::addOption("dict", "file", "Alternative blockMeshDict");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word regionName(polyMesh::defaultRegion);
    word regionPath;

    // Check if the region is specified otherwise mesh the default region
    if (args.readIfPresent("region", regionName))
    {
        Info<< nl << "Generating mesh for region " << regionName << endl;
        regionPath = regionName;
    }

    // Locate appropriate blockMeshDict
    #include "findBlockMeshDict.H"

    blockMesh blocks(meshDict, regionName);
    blocks.verbose(false);

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    const word meshInstance = runTime.constant();

    // Ensure we get information messages, even if turned off in dictionary
    blocks.verbose(true);

    autoPtr<polyMesh> meshPtr =
        blocks.mesh(IOobject(regionName, meshInstance, runTime));

    polyMesh& mesh = *meshPtr;

    // Smoothing
    const bool withQuality = args.found("write-quality");
    const bool writeSteps = args.found("writeSteps");

    label nWritten = 0;
    {
        #include "createSmoother.H"
        MeshSmoother& smoother = smootherPtr();

        if (writeSteps || withQuality)
        {
            nWritten += smoother.updateAndWrite(runTime, withQuality);
        }
        else
        {
            smoother.update();

            // Reset mesh directory to constant/polyMesh !!! Hard to find :P
            mesh.setInstance(runTime.constant());
        }
    }


    // Merge patch pairs (dictionary entry "mergePatchPairs")
    #include "mergePatchPairs.H"

    // Set any cellZones
    #include "addCellZones.H"

    // #########################################################################

    if (!nWritten)
    {
        Info<< nl << "Writing polyMesh" << endl;

        mesh.removeFiles();
        if (!mesh.write())
        {
            FatalErrorIn(args.executable())
                << "Failed writing polyMesh."
                << exit(FatalError);
        }
    }

    // #########################################################################

    //
    // write some information
    //
    {
        const polyPatchList& patches = mesh.boundaryMesh();

        Info<< "----------------" << nl
            << "Mesh Information" << nl
            << "----------------" << nl
            << "  " << "boundingBox: " << boundBox(mesh.points()) << nl
            << "  " << "nPoints: " << mesh.nPoints() << nl
            << "  " << "nCells: " << mesh.nCells() << nl
            << "  " << "nFaces: " << mesh.nFaces() << nl
            << "  " << "nInternalFaces: " << mesh.nInternalFaces() << nl;

        Info<< "----------------" << nl
            << "Patches" << nl
            << "----------------" << nl;

        forAll(patches, patchI)
        {
            const polyPatch& p = patches[patchI];

            Info<< "  " << "patch " << patchI
                << " (start: " << p.start()
                << " size: " << p.size()
                << ") name: " << p.name()
                << nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
