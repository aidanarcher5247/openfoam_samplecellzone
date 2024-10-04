/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "cellZoneSet.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "DynamicList.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "word.H"
#include "DynamicField.H"
#include "fvMesh.H"
//#include "fvcDdt.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellZoneSet, 0);
    addToRunTimeSelectionTable(sampledSet, cellZoneSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellZoneSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
    //,DynamicList<vector>& dUdtSamples // New list to hold dUdt values
) const
{
    const cellZoneMesh& zn = mesh().cellZones();
    DynamicList<label>  foundProc;
    const vectorField& my_cellC = mesh().cellCentres();
     
    //const volVectorField& U = mesh().lookupObject<volVectorField>("U"); // Access velocity field
    //volVectorField dUdt = fvc::ddt(U); // Calculate dUdt

    forAll(zn, sampleI)
    {
        if  (zn[sampleI].name() == zoneName_ )
        {   
            const cellZone& myzone  = zn[sampleI];
            
            forAll(myzone, i)
            {     
                samplingPts.append(point(my_cellC[myzone[i]]));
                samplingCells.append(myzone[i]);
                samplingFaces.append(-1);
                samplingSegments.append(0);
                samplingCurveDist.append(1.0 * i);
                
                //dUdtSamples.append(dUdt[myzone[i]]); // Add dUdt value

                foundProc.append(Pstream::myProcNo()); 
            }
                
        }
        
    }
}


void Foam::cellZoneSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;
    //DynamicList<vector> dUdtSamples; // New list for dUdt values

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
        //,dUdtSamples // Pass the new list
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();
    //dUdtSamples.shrink(); // Shrink dUdt samples

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
        //,std::move(dUdtSamples)
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::cellZoneSet::cellZoneSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    zoneName_(dict.get<word>("cellZone"))
{
    genSamples();
}


// ************************************************************************* //
