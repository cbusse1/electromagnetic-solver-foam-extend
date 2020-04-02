/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    emFoam

Description
    Solver for (quasi-)stationary electromagnetic fields produced by steady or 
    alternating currents with a given frequency. The vector potential is 
    solved in the frequency domain using a block coupled matrix solver so that 
    even high-frequency electromagnetic phenomena can be treated well. 

    Note: "Jcoil" refers to the source current producing the electromagnetic 
    field I kept this naming convention to be consistent with the documentation 
    in my master thesis, where you can read more about the theory behind the 
    solver: 
    Busse, Christian: Numerical Modeling of an Inductively Coupled 
    Plasma (ICP). Ilmenau 2019.  https://doi.org/10.22032/dbt.40314

    In the case setup the following inputs need to be specified
    in "../case/0":
        - Source current density Jcoil in A/m^2
        - Electrical conductivity sigma in A/Vm
    in "../case/constant/physicalProperties":
        - Magnetic permeability muMag (default muMag=mu0) in Vs/Am
        - Current frequency w in 1/s

    The output quantities are:
        - Magnetic vector potential A in Vs/m
        - Magnetic flux density B in Vs/m^2
        - Magnetic field strength H in A/m
        - Induced current density Jind in A/m^2
        - (Time-averaged) Joule heat density in W/m^3  
        - (Time-averaged) Lorentz-force density fL in N/m^3

Author
    Christian Busse

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H" // load block coupled matrix solver

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initConvergenceCheck.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
/*
#       include "readBlockSolverControls.H" // I didn't use these, but someone might want to have it
#       include "readFieldBounds.H"
*/
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update auxiliary variable DA 
        DA = muMag*sigma*w;

        // Initialize the vector potential matrix
        fvBlockMatrix<vector6> AEqn(A);

            fvVectorMatrix AREqn
            (
                fvm::laplacian(AR) == - muMag*Jcoil
            ); 

            fvVectorMatrix AIEqn
            (
                 fvm::laplacian(AI) 
            );

        //insert fvVectorMatrix equations into the fvBlockMatrix
        AEqn.insertEquation(0, AREqn);
        AEqn.insertEquation(3, AIEqn);

        // Add off-diagonal coupling terms
        AEqn.insertEquationCoupling(0, 3, DA);
        AEqn.insertEquationCoupling(1, 4, DA);
        AEqn.insertEquationCoupling(2, 5, DA);
        AEqn.insertEquationCoupling(3, 0, -DA);
        AEqn.insertEquationCoupling(4, 1, -DA);
        AEqn.insertEquationCoupling(5, 2, -DA);

        // Solve the block matrix
        maxResidual = cmptMax(AEqn.solve().initialResidual());

        // Retrieve solution of the vector potential A
        AEqn.retrieveSolution(0, AR.internalField());
        AEqn.retrieveSolution(3, AI.internalField());

        // Calculate the vector potential magnitude
        Amag = sqrt((AR & AR) + (AI & AI)); // & - dot product of two vectors

        // Retrieve B and H fields
        BR = fvc::curl(AR);
        BI = fvc::curl(AI);
        Bmag = sqrt((BR & BR) + (BI & BI)); // & - dot product of two vectors        
        Hmag = Bmag/muMag;   

        // Resolve the induced current density J in [A/m^2]
        JR = sigma*w*AI;  
        JI = -sigma*w*AR;
        Jind = sqrt((JR & JR) + (JI & JI)); // & - dot product of two vectors  
        
        // Compute the time-averaged Joule heat (power density) and Lorentz-force 
        qJ = sigma/2.0 * sqr(w) * sqr(Amag); // Dissipated power density in [W/m^3]   
        fL = 0.5 * ((JR ^ BR) + (JI ^ BI)); // Lorentz-force in [N/m^3]

        // Compute the total dissipated power in [W]
        dimensionedScalar QJ = fvc::domainIntegrate(qJ);       

        Info<< "----- Total dissipated power: " << nl << endl;
        Info<< "QJ  = " << QJ.value() << nl << endl;

        Info<< "EM Solver: ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        #include "convergenceCheck.H"

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
