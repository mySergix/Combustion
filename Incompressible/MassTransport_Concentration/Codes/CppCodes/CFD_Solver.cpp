//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS                                      //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/CFD_Solver.h"

CFD_Solver::CFD_Solver(Memory M1, ReadData R1, Parallel P1){

    //Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	NX = R1.ProblemNumericalData[2];
	NY = R1.ProblemNumericalData[3];
	NZ = R1.ProblemNumericalData[4];
	Halo = 2;
	HP = 2;

    //Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	Dh = (4.0 * Zdominio * Ydominio) / (2.0 * (Zdominio + Ydominio));
	gamma_geometry = Zdominio / Ydominio;

	//Datos necesarios para computación paralela
	Rango = P1.Rango;
	Procesos = P1.Procesos;
	Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

    for (int i = 0; i < Procesos; i++){
        Ix[i] = P1.Ix[i];
        Fx[i] = P1.Fx[i];
    }

    //Datos Físicos del Problema
	Rho = R1.ProblemPhysicalData[0];
	Uref = R1.ProblemPhysicalData[1];
	Reynolds = R1.ProblemPhysicalData[2];
	
	Rayleigh = R1.ProblemPhysicalData[3];
	Cp = R1.ProblemPhysicalData[4];
	Prandtl = R1.ProblemPhysicalData[5];

	U.Gravity = R1.ProblemPhysicalData[6];
	V.Gravity = R1.ProblemPhysicalData[7];
	W.Gravity = R1.ProblemPhysicalData[8];

	//mu = (Uref * Xdominio) / Reynolds;
	
	Tleft = 1.0;
	Tright = 0.0;

	To = 0.50 * (Tleft + Tright);

	Beta = 1.0/To;	
	Producto = (pow(Rho,2)*abs(V.Gravity)*pow(Xdominio,3)*Beta*abs(Tleft-Tright)*Prandtl)/Rayleigh;

	if(Problema == 2){
		mu = sqrt(Producto);
	}
	K = (Cp*mu)/Prandtl;
    
	
	//nu = 15.24e-6; // m2/s
    //mu = nu * Rho;
	//(1.0 * Dh) / Reynolds;
	//;
/*
	Schmidt = 0.59;
	D_AB = nu / Schmidt;
	w_av = (mu * Reynolds) / (Rho * Dh);
	K = (Cp*mu)/Prandtl;
	Co = 0.531;
	MW_H2O = 18.0;
	hfg = (40.7 * 1000.0) / (1e-3 * MW_H2O);
	
	T_water = 16.0; // Celsius
	To = 22.4; // Celsius

	Beta = 0.00343;
	Beta_Concentration = 0.513;
	*/
    ConvergenciaGS = R1.ProblemData[3];
	ConvergenciaGlobal = R1.ProblemData[4];

}

// Files of the class (CFD)
#include "Matrix_Index.cpp"
#include "CFD_Solver_Memory.cpp"
#include "CFD_Solver_Utilities.cpp"
#include "CFD_Solver_BoundaryConditions_Differentially.cpp"
#include "CFD_Solver_PoissonCoeffs.cpp"
#include "CFD_Solver_MomentumDiffusion.cpp"
#include "CFD_Solver_MomentumConvection.cpp"
#include "CFD_Solver_Energy.cpp"
#include "CFD_Solver_RunSolver.cpp"

// Files of the class (Species)
#include "Species_Solver_Memory.cpp"
#include "Species_Solver_BoundaryConditions.cpp"
#include "Species_Solver_Diffusion.cpp"
#include "Species_Solver_Convection.cpp"