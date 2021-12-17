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

// Function to calculate the binary diffusion based on Champan-Enskog model
double Species_Solver::Get_BinaryDiff_ChampanEnskog(int SP, double Temperature, double Pressure, int i, int j, int k){
int i;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn;

    for (i = 0; i < N_Species; i++){
        if (i =! SP){
            W_ab = power(1.0 / Species[SP].W_molar + 1.0 / Species[i].W_molar, - 1.0);
            Sigma_ab = (Species[SP].sigma + Species[i].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[SP].Epsilon * Species[i].Epsilon);
            E_ab = Epsilon_ab / kB;
            Tn = Temperature / E_ab;
            OmegaD = 1.06036 / power(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_ab = 10.1325 * (0.001858 * power(Tn, 1.5) * power(W_ab, -0.5)) / (Pressure * power(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[i].X[LM(i,j,k,0)] / D_ab;
        }    
    }
       
    return Sum;

}

// Function to calculate the binary difussion based on Wilke-Lee model 
double Species_Solver::Get_BinaryDiff_WilkeLee(int SP, double Temperature, double Pressure, int i, int j, int k){
int i;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn;

    for (i = 0; i < N_Species; i++){
        if (i =! SP){
            W_ab = power(1.0 / Species[SP].W_molar + 1.0 / Species[i].W_molar, - 1.0);
            Sigma_ab = (Species[SP].sigma + Species[i].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[SP].Epsilon * Species[i].Epsilon);
            E_ab = Epsilon_ab / kB;
            Tn = Temperature / E_ab;
            OmegaD = 1.06036 / power(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_ab = 10.1325 * ((0.0027 - 0.0005 * power(W_ab, - 0.5)) * power(Tn, 1.5) * power(W_ab, -0.5)) / (Pressure * power(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[i].X[LM(i,j,k,0)] / D_ab;
        }    
    }
       
    return Sum;

}

// Function to calculate the diffusion coefficient of a species
void Species_Solver::Get_DiffusionCoefficient_FickModel(int SP, CFD_Solver CFD_S1){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].D_am[LM(i,j,k,0)] = (1.0 - Species[SP].X[LM(i,j,k,0)]) / Get_BinaryDif_ChampanEnskog(SP, CFD_S1.T.Pres[LP(i,j,k,0)], CFD_S1.Pressure.Pres[LP(i,j,k,0)], i, j, k);
            }
        }
    }

}

// Function to calculate the walls value of diffusion coefficients
void Species_Solver::Get_WallsDiffusionCoefficients(int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].D_am_Wall_U[LMU(i,j,k,0)] = 0.50 * (Species[SP].D_am[LM(i-1,j,k,0)] + Species[SP].D_am[LM(i,j,k,0)]);
                // Faltan direcciones Y y Z
            }
        }
    }

}

// Function to calculate the diffusion term of a species
void Species_Solver::Get_SpeciesDiffusion(Mesher MESH, CFD_Solver CFD_S1, int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].DiffusionTerm[LM(i,j,k,0)] = (1.0 / MESH.VolMP) * (
                                                       - MESH.Surf[LM(i,j,k,0)] * CFD_S1.Density.Wall_U[LMU(i,j,k,0)] * Species[SP].D_am_Wall_U[LMU(i,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i-1,j,k,0)] + MESH.DeltaP[LM(i,j,k,0)])) * (Species[SP].Y[LM(i,j,k,0)] - Species[SP].Y[LM(i-1,j,k,0)])
                                                       + MESH.Surf[LM(i,j,k,0)] * CFD_S1.Density.Wall_U[LMU(i+1,j,k,0)] * Species[SP].D_am_Wall_U[LMU(i+1,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k,0)] + MESH.DeltaP[LM(i+1,j,k,0)])) * (Species[SP].Y[LM(i+1,j,k,0)] - Species[SP].Y[LM(i,j,k,0)])
                                                       - MESH.Surf[LM(i,j,k,1)] * CFD_S1.Density.Wall_V[LMV(i,j,k,0)] * Species[SP].D_am_Wall_V[LMV(i,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j-1,k,1)] + MESH.DeltaP[LM(i,j,k,1)])) * (Species[SP].Y[LM(i,j,k,0)] - Species[SP].Y[LM(i,j-1,k,0)])
                                                       + MESH.Surf[LM(i,j,k,1)] * CFD_S1.Density.Wall_V[LMV(i,j+1,k,0)] * Species[SP].D_am_Wall_V[LMV(i,j+1,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k,1)] + MESH.DeltaP[LM(i,j+1,k,1)])) * (Species[SP].Y[LM(i,j+1,k,0)] - Species[SP].Y[LM(i,j,k,0)])
                                                       - MESH.Surf[LM(i,j,k,2)] * CFD_S1.Density.Wall_W[LMW(i,j,k,0)] * Species[SP].D_am_Wall_W[LMW(i,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k-1,2)] + MESH.DeltaP[LM(i,j,k,2)])) * (Species[SP].Y[LM(i,j,k,0)] - Species[SP].Y[LM(i,j,k-1,0)])
                                                       + MESH.Surf[LM(i,j,k,2)] * CFD_S1.Density.Wall_W[LMW(i,j,k+1,0)] * Species[SP].D_am_Wall_W[LMW(i,j,k+1,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k,2)] + MESH.DeltaP[LM(i,j,k+1,2)])) * (Species[SP].Y[LM(i,j,k+1,0)] - Species[SP].Y[LM(i,j,k,0)])
                                                       );
            }
        }
    }

}

// Function to calculate the convection term of a species
void Species_Solver::Get_SpeciesConvection(Mesher MESH, CFD_Solver CFD_S1, int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].ConvectiveTerm[LM(i,j,k,0)] = (1.0/MESH.Vol[LM(i,j,k,0)]) * (
                                                     + MESH.Surf[LM(i,j,k,0)] * CFD_S1.U.Wall_U[LMU(i+1,j,k,0)] * Species[SP].Y_Wall_U[LMU(i+1,j,k,0)]
                                                     - MESH.Surf[LM(i,j,k,0)] * CFD_S1.U.Wall_U[LMU(i,j,k,0)] * Species[SP].Y_Wall_U[LMU(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,1)] * CFD_S1.V.Wall_V[LMV(i,j+1,k,0)] * Species[SP].Y_Wall_V[LMV(i,j+1,k,0)]
                                                     - MESH.Surf[LM(i,j,k,1)] * CFD_S1.V.Wall_V[LMV(i,j,k,0)] * Species[SP].Y_Wall_V[LMV(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,2)] * CFD_S1.W.Wall_W[LMW(i,j,k+1,0)] * Species[SP].Y_Wall_W[LMW(i,j,k+1,0)]
                                                     - MESH.Surf[LM(i,j,k,2)] * CFD_S1.W.Wall_W[LMW(i,j,k,0)] * Species[SP].Y_Wall_W[LMW(i,j,k,0)]
                                                     );
            }
        }
    }

}

// Function to calculate the step contribution to each species
void Species_Solver::Get_StepContribution_Species(CFD_Solver CFD_S1, int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].ContributionPres[LM(i,j,k,0)] = (1.0 / CFD_S1.Density.Pres[LM(i,j,k,0)]) * (- Species[SP].ConvectiveTerm[LM(i,j,k,0)] - Species[SP].DiffusionTerm[LM(i,j,k,0)]);
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
                                              + DeltaT * ((1.0 + Beta) * Species[SP].ContributionPres[LM(i,j,k,0)] - Beta * Species[SP].ContributionPast[LM(i,j,k,0)]); 
            }
        }
    }

}

// Function to update the species amp fields
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
                    Sum += Species[sp].Y_Pres[LM(i,j,k,0)] / Species[sp].W_molar;
                }
                W_molar_Average = 1.0 / Sum;    

                for (sp = 0; sp < N_Species; sp++){
                    Species[sp].X[LM(i,j,k,0)] = (W_molar_Average / Species[sp].W_molar) * Species[sp].Y_Pres[LM(i,j,k,0)];
                }

            }
        }
    }
}