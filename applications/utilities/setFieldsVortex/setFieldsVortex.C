/*---------------------------------------------------------------------------*\
This file was modified or created at the National Renewable Energy
Laboratory (NREL) on January 6, 2012 in creating the SOWFA (Simulator for
Offshore Wind Farm Applications) package of wind plant modeling tools that
are based on the OpenFOAM software. Access to and use of SOWFA imposes
obligations on the user, as set forth in the NWTC Design Codes DATA USE
DISCLAIMER AGREEMENT that can be found at
<http://wind.nrel.gov/designcodes/disclaimer.html>.
\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    setFieldsVortex

Description
    Initializes the flow field for a Lamb-Oseen Vortex.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "interpolateXY.H"
#include "interpolateSplineXY.H"
#include "Random.H"
#include "wallDist.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Read in the existing solution files.   
Info << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Get access to the input dictionary.
const word dictName("setFieldsVortexDict");
#include "setSystemMeshDictionaryIO.H"
Info << "Reading " << dictPath << endl;
IOdictionary dict(dictIO);
/*
IOdictionary dict
(
    IOobject
    (
        "dict",
        runTime.time().system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
*/



// Read in the dict entries.
scalar xMin(dict.lookupOrDefault<scalar>("xMin",0.0));
scalar yMin(dict.lookupOrDefault<scalar>("yMin",0.0));
bool updateInternalFields(dict.lookupOrDefault<bool>("updateInternalFields",true));
bool updateBoundaryFields(dict.lookupOrDefault<bool>("updateBoundaryFields",false));

// New constants for the Lamb-Oseen Vortex
scalar xc(dict.lookupOrDefault<scalar>("xc",2500.0));
scalar yc(dict.lookupOrDefault<scalar>("yc",2500.0));
scalar Uo(dict.lookupOrDefault<scalar>("Uo",50.0));
scalar Rpeak(dict.lookupOrDefault<scalar>("Rpeak",1000.0));

// Read the convecting wind speed
scalar Uc(dict.lookupOrDefault<scalar>("Uc",0.0));
scalar Vc(dict.lookupOrDefault<scalar>("Vc",0.0));


// Update the interior fields.
if (updateInternalFields)
{
    // Velocity.
    Info << "Updating internal U field..." << endl;
    forAll(U,cellI)
    {
        scalar x = mesh.C()[cellI].x() - xMin;
        scalar y = mesh.C()[cellI].y() - yMin;
        
	// Calculate the local radius and tangential velocity
	// Tangential velocity includes a peak velocity radius correction from Davenport et al. 1996 (12% offset)
	scalar r = Foam::pow(Foam::pow((x - xc),2) + Foam::pow((y - yc),2),0.5);
	scalar u_theta = Uo*(1+1/(2*1.26))*(Rpeak/r)*(1 - Foam::exp(-1.26*Foam::pow((r/Rpeak),2)));
	
	// Calculate the local theta
	scalar thetaU = Foam::asin((y - yc)/r);
	scalar thetaV = Foam::acos((x - xc)/r);

       	// Compute velocity components for the Lamb-Oseen Vortex
        U[cellI].x() = -u_theta*Foam::sin(thetaU) + Uc;
	U[cellI].y() =  u_theta*Foam::cos(thetaV) + Vc;
     }
}


// Update the boundary field.
if (updateBoundaryFields)
{
    Info << "Updating boundaries..." << endl;
    U.correctBoundaryConditions();
}


// Write out the updated fields.
Info<< "Writing field U" << endl;
U.write();


return 0;
}

// ************************************************************************* //

