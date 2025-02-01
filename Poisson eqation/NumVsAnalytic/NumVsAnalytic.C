/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    NumVsAnalytic

Description
    post-process of comparison of analytical and numerical solutions.

\*---------------------------------------------------------------------------*/


#include "NumVsAnalytic.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace functionObjects
  {
      defineTypeNameAndDebug(NumVsAnalytic, 0);
      addToRunTimeSelectionTable(functionObject, NumVsAnalytic, dictionary);
  }
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * // 

Foam::functionObjects::NumVsAnalytic::NumVsAnalytic
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    time_(runTime),
    mesh_(runTime.lookupObject<fvMesh>(polyMesh::defaultRegion)),
	
    phi_
    (
        IOobject
		(
          "phi",
          runTime.timeName(),
          mesh_,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
		),
        mesh_,
        dimensionedScalar 
       (
          "phi",     
          dimensionSet(0, 2, 0, 0, 0, 0, 0), 
          0.0       
       )
   ),
   phiAnalytic_
    (
      IOobject
      (
          "phiAnalytic",
          runTime.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar 
      (
          "phiAnalytic",     
          dimensionSet(0, 2, 0, 0, 0, 0, 0), 
          0.0       
      )
   ),
   err_
   (
      IOobject
      (
          "err",
          runTime.timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mesh_,
      dimensionedScalar 
      (
          "err",     
          dimensionSet(0, 2, 0, 0, 0, 0, 0), 
          0.0       
      )
   )
{
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //
bool Foam::functionObjects::NumVsAnalytic::execute()
{
    const scalar t = time_.time().value(); 
    
	forAll(phiAnalytic_, cell)
	{
		const scalar x = mesh_.C()[cell].x();
		const scalar y = mesh_.C()[cell].y();
		
		phiAnalytic_.ref()[cell] = Foam::cos(x * t) * Foam::cos(y * t) + 1.5 * Foam::sin(t) * Foam::sin(t) * x * x * y * y;
	}
	err_.ref() = phi_ - phiAnalytic_;
	
	report();
	
    return true;
}

void Foam::functionObjects::NumVsAnalytic::report()
{
	const scalar aver_err = Foam::sqrt(gSum(mesh_.V() * magSqr(err_)) / gSum(mesh_.V()));
	
    Info << "Average error = " << aver_err << endl; 
}

bool Foam::functionObjects::NumVsAnalytic::start()
{ 
    return execute();
}

bool Foam::functionObjects::NumVsAnalytic::end()
{
    return execute();
}

void Foam::functionObjects::NumVsAnalytic::updateMesh(const mapPolyMesh& map)
{
    execute(); 
}

void Foam::functionObjects::NumVsAnalytic::movePoints(const polyMesh& mesh)
{
    execute();  
}

bool Foam::functionObjects::NumVsAnalytic::write()
{
    if (time_.writeTime())
    {
        phiAnalytic_.write();
        err_.write();
    }
    return true;
}
