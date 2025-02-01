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
    PoissonSolver

Description
    Poisson equation solver.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Poisson equation solver."
    );

	  #include "dimensionedScalar.H" 
	  #include "volFields.H"        
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
	
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info << endl << "Time = " << runTime.timeName() << endl;

        // Создание правой части
        const scalar t = runTime.time().value(); 
        scalar fSum = 0.0;
		   
        forAll(f.ref(), cell)
        {
            const scalar x = mesh.C()[cell].x();
            const scalar y = mesh.C()[cell].y();
            
            f.ref()[cell] = - 2.0 * t * t * Foam::cos(x * t) * Foam::cos(y * t) + 3 * Foam::sin(t) * Foam::sin(t) * (x * x + y * y);
            fSum += f.ref()[cell] * mesh.V()[cell];
        }
        
        // Коррекция f для выполнения ∫fdV = 0
        scalar integr_val = fSum / gSum(mesh.V());
        f.ref() -= integr_val;  

        // Решение уравнения Пуассона
        fvScalarMatrix phiEqn(fvm::laplacian(phi) == f);
        phiEqn.solve();		
        
        runTime.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
