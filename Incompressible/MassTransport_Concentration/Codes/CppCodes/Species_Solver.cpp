//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SPECIES SOLVER CLASS                                  //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/CFD_Solver.h"
#include "../HeaderCodes/Species_Solver.h"

Species_Solver::Species_Solver(Memory M1, ReadData R1, Parallel P1){

    // Data
    NX = R1.ProblemNumericalData[2];
	NY = R1.ProblemNumericalData[3];
	NZ = R1.ProblemNumericalData[4];
	Halo = 2;
	HP = 2;

    // Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	//Datos necesarios para computación paralela
	Rango = P1.Rango;
	Procesos = P1.Procesos;
	Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

    for (int i = 0; i < Procesos; i++){
        Ix[i] = P1.Ix[i];
        Fx[i] = P1.Fx[i];
    }

	nu = 15.24e-6; // m2/s
    mu = nu * Rho;

	Schmidt = 0.59;
	D_AB = nu / Schmidt;

    Co = 0.531;
	MW_H2O = 18.0;
	hfg = (40.7 * 1000.0) / (1e-3 * MW_H2O);
	T_water = 16.0;

}

// Files of the class
#include "Matrix_Index.cpp"
#include "Species_Solver_Memory.cpp"
#include "Species_Solver_BoundaryConditions.cpp"
#include "Species_Solver_Diffusion.cpp"
#include "Species_Solver_Convection.cpp"

// Function to calculate the concentration contribution
void Species_Solver::Get_Concentration(){
int i, j, k;

    // Species Contributions
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				Species[0].ContributionPres[LP(i,j,k,0)] = Species[0].Diffusive[LP(i,j,k,0)] - Species[0].Convective[LP(i,j,k,0)];
			}
		}
	}

    // Temperature T Calculations
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				Species[0].C_Fut[LP(i,j,k,0)] = Species[0].C_Pres[LP(i,j,k,0)] + DeltaT*(1.50*Species[0].ContributionPres[LP(i,j,k,0)] - 0.50*Species[0].ContributionPast[LP(i,j,k,0)]);
			}
		}
	}

}