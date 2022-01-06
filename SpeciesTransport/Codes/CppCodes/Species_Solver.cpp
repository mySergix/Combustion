//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SPECIES SOLVER CLASS                                  //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcess.h"
#include "../HeaderCodes/CFD_Solver.h"
#include "../HeaderCodes/Species_Solver.h"



using namespace std;

#define LM(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo) + ((Fx[Rango] - Ix[Rango] + 2*Halo) * (NY + 2*Halo) * (NZ + 2*Halo)) * (dim)
#define GM(i,j,k,dim) (NY + 2*Halo) * (NZ + 2*Halo) * ((i) + Halo) + (NZ + 2*Halo) * ((j) + Halo) + ((k) + Halo) + (NY + 2*Halo) * (NZ + 2*Halo) * (NX + 2*Halo) * (dim)

Species_Solver::Species_Solver(Memory M1, ReadData R1, Parallel P1){
	
    // Data necessary

        // Geometry data
        Xdominio = R1.GeometryData[0];
        Ydominio = R1.GeometryData[1];
        Zdominio = R1.GeometryData[2];

        // Meshing data
        NX = R1.ProblemNumericalData[0];
        NY = R1.ProblemNumericalData[1];
        NZ = R1.ProblemNumericalData[2];

        // Parallel computing data
        Rango = P1.Rango;
        Procesos = P1.Procesos;

        Ix = M1.AllocateInt(Procesos, 1, 1, 1);
        Fx = M1.AllocateInt(Procesos, 1, 1, 1);

        for (int i = 0; i < Procesos; i++){
            Ix[i] = P1.Ix[i];
            Fx[i] = P1.Fx[i];
        }

        Halo = 2;

        kB = 1.380649e-23; // Boltzman constant

}

// Parts of the Species Solver Class
#include "Species_Solver_Memory.cpp"
#include "Species_Solver_JANAF.cpp"
#include "Species_Solver_Diffusion.cpp"

// Function to calculate the convective term of each species
void Species_Solver::Get_SpeciesConvection(Mesher MESH, CFD_Solver CFD_S1, int SP){
int i, j, k;

    // Importante, la aproximacion solo vale para flujos sin fuerzas volumetricas y SIN GRADIENTES DE PRESION

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].ContributionPres[LM(i,j,k,0)] = (1.0 / MESH.VolMP) * (
                                                       - MESH.Surf[LM(i,j,k,0)] * (CFD_S1.Density.Wall_U[LMU(i,j,k,0)] + Species[SP].U_Diff[LMU(i,j,k,0)]) * Species[SP].Y_Wall_U[LMU(i,j,k,0)]
                                                       + MESH.Surf[LM(i,j,k,0)] * (CFD_S1.Density.Wall_U[LMU(i+1,j,k,0)] + Species[SP].U_Diff[LMU(i+1,j,k,0)]) * Species[SP].Y_Wall_U[LMU(i+1,j,k,0)]
                                                       - MESH.Surf[LM(i,j,k,1)] * (CFD_S1.Density.Wall_V[LMV(i,j,k,0)] + Species[SP].V_Diff[LMV(i,j,k,0)]) * Species[SP].Y_Wall_V[LMV(i,j,k,0)] 
                                                       + MESH.Surf[LM(i,j,k,1)] * (CFD_S1.Density.Wall_V[LMV(i,j+1,k,0)] + Species[SP].V_Diff[LMV(i,j+1,k,0)]) * Species[SP].Y_Wall_V[LMV(i,j+1,k,0)]
                                                       - MESH.Surf[LM(i,j,k,2)] * (CFD_S1.Density.Wall_W[LMW(i,j,k,0)] + Species[SP].W_Diff[LMW(i,j,k,0)]) * Species[SP].Y_Wall_W[LMW(i,j,k,0)]
                                                       + MESH.Surf[LM(i,j,k,2)] * (CFD_S1.Density.Wall_W[LMW(i,j,k+1,0)] + Species[SP].W_Diff[LMW(i,j,k+1,0)]) * Species[SP].Y_Wall_W[LMW(i,j,k+1,0)]
                                                       );
            }
        }
    }

}

// Function to integrate the species equation
void Species_Solver::Get_TemporalIntegration_Species(int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].Y_Fut[LM(i,j,k,0)] = (2.0 * Beta * Species[SP].Y_Pres[LM(i,j,k,0)] - (Beta - 0.50) * Species[SP].Y_Past[LM(i,j,k,0)]) / (Beta + 0.50) 
                                              + DeltaT * ((1.0 + Beta) * (1.0 / CFD_S1.Density.Pres[LM(i,j,k,0)]) * Species[SP].ContributionPres[LM(i,j,k,0)] - Beta * (1.0 / CFD_S1.Density.Past[LM(i,j,k,0)]) * Species[SP].ContributionPast[LM(i,j,k,0)]); 
            }
        }
    }

}

// Function to update the species properties fields
void Species_Solver::Get_Update(int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].Y_Past[LM(i,j,k,0)] = Species[SP].Y_Pres[LM(i,j,k,0)];
                Species[SP].Y_Pres[LM(i,j,k,0)] = Species[SP].Y_Fut[LM(i,j,k,0)];

                Species[SP].ContributionPast[LM(i,j,k,0)] = Species[SP].ContributionPres[LM(i,j,k,0)];
            }
        }
    }

}

// Function to calculate the molar fraction of each species
void Species_Solver::Get_MolarFraction_X(){
int i, j, k, sp;
double Sum;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){

                Sum = 0.0;
                for (sp = 0; sp < N_Species; sp++){
                    Sum += Species[sp].Y_Pres[LM(i,j,k,0)] / Species[sp].Wmolar;
                }
                W_molar_Average = 1.0 / Sum;    

                for (sp = 0; sp < N_Species; sp++){
                    Species[sp].X[LM(i,j,k,0)] = (W_molar_Average / Species[sp].Wmolar) * Species[sp].Y_Pres[LM(i,j,k,0)];
                }

            }
        }
    }
    
}