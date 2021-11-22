//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS                                      //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/CFD_Solver.h"
#include "../HeaderCodes/Species_Solver.h"
#include "../HeaderCodes/PostProcess.h"


using namespace std;

#define LM(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo) + ((Fx[Rango] - Ix[Rango] + 2*Halo) * (NY + 2*Halo) * (NZ + 2*Halo)) * (dim)
#define LMU(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo)
#define LMV(i,j,k,dim) ((NY + 1 + 2*Halo) * (NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo)
#define LMW(i,j,k,dim) ((NY + 2*Halo) * (NZ + 1 + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 1 + 2*Halo) * ((j) + Halo)) + ((k) + Halo)

#define GM(i,j,k,dim) (NY + 2*Halo) * (NZ + 2*Halo) * ((i) + Halo) + (NZ + 2*Halo) * ((j) + Halo) + ((k) + Halo) + (NY + 2*Halo) * (NZ + 2*Halo) * (NX + 2*Halo) * (dim)

#define TOP_BOT(i,j,k,dim) ((NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((k) + Halo)
#define LEF_RIG(i,j,k,dim) ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo)
#define HER_THE(i,j,k,dim) ((NY + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((j) + Halo)

CFD_Solver::CFD_Solver(Memory M1, ReadData R1, Parallel P1){

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

};

// Function to allocate memory for all matrix
void CFD_Solver::AllocateMemory(Memory M1){

    // Mass Conservation
    Density = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Rho_W.Wall_U = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1); 
    Rho_W.Wall_V = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Rho_W.Wall_W = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);
    Density_Convective = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    Rho_Boun.Bottom = M1.AllocateDouble(Fx[Rango] - Ix[Rango], 1, NZ, 1);
    Rho_Boun.Top = M1.AllocateDouble(Fx[Rango] - Ix[Rango], 1, NZ, 1);

    Rho_Boun.Here = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, 1, 1);
    Rho_Boun.There = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, 1, 1);

    if (Rango == 0){
        Rho_Boun.Left = M1.AllocateDouble(Fx[Rango] - Ix[Rango], 1, NZ, 1);
        
    }
    if (Rango == Procesos - 1){
        Rho_Boun.Right = M1.AllocateDouble(Fx[Rango] - Ix[Rango], 1, NZ, 1);
    }
    
}

// Convective Scheme Function
double CFD_Solver::CS(double X, double V, double X1, double Phi1, double X2, double Phi2, double X3, double Phi3, double X4, double Phi4){
double XD, PHID, XC, PHIC, XU, PHIU;
double PhiC_, XC_, XE_, PHIE_, PhiE;

    if (V >= 0){
        XD = X3;
        PHID = Phi3;
        XC = X2;
        PHIC = Phi2;
        XU = X1;
        PHIU = Phi1;
    }
    else{
        XD = X2;
        PHID = Phi2;
        XC = X3;
        PHIC = Phi3;
        XU = X4;
        PHIU = Phi4;
    }

    // Non - Dimensionalization
    PhiC_ = (PHIC - PHIU) / (PHID - PHIU);
    XC_ = (XC - XU) / (XD - XU);
    XE_ = (X - XU) / (XD - XU);

    // CDS Scheme
    PHIE_ = ((XE_ - XC_) / (1 - XC_)) + ((XE_ - 1) / (XC_ - 1)) * PhiC_;

    // Dimensionalization
    PhiE = PHIU + (PHID - PHIU) * PHIE_;

    return PhiE;

}

// Function to get the boundary conditions of each property
void CFD_Solver::Get_BoundaryConditions(){
int i, j, k;

    // Parte Top y Bottom
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (k = - Halo; k < NZ + Halo; k++){
            // Parte Bottom
            Rho_Boun.Bottom[TOP_BOT(i,0,k,0)] = 0.0;

            //Parte Top
            Rho_Boun.Top[TOP_BOT(i,NY,k,0)] = 0.0;
        }
    }

    // Parte Here y There
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            // Parte Here
            Rho_Boun.Here[HER_THE(i,j,0,0)] = 0.0;

            // Parte There
            Rho_Boun.There[HER_THE(i,j,NZ,0)] = 0.0;
        }
    }

    // Parte Left
    if (Rango == 0){

        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Rho_Boun.Left[LEF_RIG(0,j,k,0)] = 1.0;
            }
        }

    }

    // Parte Right
    if (Rango == Procesos - 1){

        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Rho_Boun.Left[LEF_RIG(0,j,k,0)] = Density[LM(NX-1,j,k,0)];
            }
        }
        
    }

}

// Function to communicate all the velocity fields (local fields)
void CFD_Solver::CommunicateVelocities(Parallel P1, double *UFIELD, double *VFIELD, double *WFIELD){

    P1.CommunicateLocalMatrix(UFIELD, UFIELD); // Velocity U
    P1.CommunicateLocalMatrix(VFIELD, VFIELD); // Velocity V
    P1.CommunicateLocalMatrix(WFIELD, WFIELD); // Velocity W

}

// Calculation of Property Value at the Walls of the Volumes
void CFD_Solver::Get_WallsProperty_Scalar(Parallel P1, double *Mesh, double *PropertyField, double *UFIELD, double *VFIELD, double *WFIELD, Walls &WallProperty, Boundaries &Bound){
int i, j, k;

    // Communication of the local metrix of the properties
    P1.CommunicateLocalMatrix(PropertyField, PropertyField);

    // MU Walls
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                WallProperty.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (Mesh[LM(i - 1,j,k,0)] + Mesh[LM(i,j,k,0)]), 0.50 * (UFIELD[LM(i - 1,j,k,0)] + UFIELD[LM(i,j,k,0)]), Mesh[LM(i-2,j,k,0)], PropertyField[LM(i-2,j,k,0)], Mesh[LM(i-1,j,k,0)], PropertyField[LM(i-1,j,k,0)], Mesh[LM(i,j,k,0)], PropertyField[LM(i,j,k,0)], Mesh[LM(i+1,j,k,0)], PropertyField[LM(i+1,j,k,0)]);
            }
        }
    }

    // MV Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                WallProperty.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j-1,k,1)] + Mesh[LM(i,j,k,1)]), 0.50 * (VFIELD[LM(i,j-1,k,0)] + VFIELD[LM(i,j,k,0)]), Mesh[LM(i,j-2,k,1)], PropertyField[LM(i,j-2,k,0)], Mesh[LM(i,j-1,k,1)], PropertyField[LM(i,j-1,k,0)], Mesh[LM(i,j,k,1)], PropertyField[LM(i,j,k,0)], Mesh[LM(i,j+1,k,1)], PropertyField[LM(i,j+1,k,0)]);
            }
        }
    }

    // MW Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 1; k < NZ; k++){
                WallProperty.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j,k-1,2)] + Mesh[LM(i,j,k,2)]), 0.50 * (WFIELD[LM(i,j,k-1,0)] + WFIELD[LM(i,j,k,0)]), Mesh[LM(i,j,k-2,2)], PropertyField[LM(i,j,k-2,0)], Mesh[LM(i,j,k-2,2)], PropertyField[LM(i,j-1,k,0)], Mesh[LM(i,j,k,2)], PropertyField[LM(i,j,k,0)], Mesh[LM(i,j,k+1,2)], PropertyField[LM(i,j,k+1,0)]);
            }
        }
    }

    ApplyBoundaries(WallProperty, Bound);   

}

// Function for the application of the boundary condition to the walls of the volumes
void CFD_Solver::ApplyBoundaries(Walls &WallProperty, Boundaries &Bound){
int i, j, k;

    // Left Side
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){               
                WallProperty.Wall_U[LMU(0,j,k,0)] = Bound.Left[LEF_RIG(0,j,k,0)];
            }
        }
    }

    // Right Side
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){               
                WallProperty.Wall_U[LMU(NX,j,k,0)] = Bound.Left[LEF_RIG(NX,j,k,0)];
            }
        }
    }

    // Top and Bottom Sides
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            // Bottom Part
            WallProperty.Wall_V[LMV(i,0,k,0)] = Bound.Bottom[TOP_BOT(i,0,k,0)];

            // Top Part
            WallProperty.Wall_V[LMV(i,NY,k,0)] = Bound.Top[TOP_BOT(i,NY,k,0)];
        }
    }

    // Here and There Sides
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            // Here Part
            WallProperty.Wall_W[LMW(i,j,0,0)] = Bound.Here[HER_THE(i,j,0,0)];

            // There Part
            WallProperty.Wall_W[LMW(i,j,NZ,0)] = Bound.There[HER_THE(i,j,NZ,0)];
        }
    }

}

// Calculation of Velocities at the Walls of the Volumes
void CFD_Solver::Get_WallsVelocities(double *Mesh, double *UFIELD, double *VFIELD, double *WFIELD, Walls &WallProperty1, Boundaries &Bound1, Walls &WallProperty2, Boundaries &Bound2, Walls &WallProperty3, Boundaries &Bound3){
int i, j, k;

    // MU Walls
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                WallProperty1.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (Mesh[LM(i - 1,j,k,0)] + Mesh[LM(i,j,k,0)]), 0.50 * (UFIELD[LM(i - 1,j,k,0)] + UFIELD[LM(i,j,k,0)]), Mesh[LM(i-2,j,k,0)], UFIELD[LM(i-2,j,k,0)], Mesh[LM(i-1,j,k,0)], UFIELD[LM(i-1,j,k,0)], Mesh[LM(i,j,k,0)], UFIELD[LM(i,j,k,0)], Mesh[LM(i+1,j,k,0)], UFIELD[LM(i+1,j,k,0)]);
            }
        }
    }

    ApplyBoundaries(WallProperty1, Bound1);  

    // MV Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                WallProperty2.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j-1,k,1)] + Mesh[LM(i,j,k,1)]), 0.50 * (VFIELD[LM(i,j-1,k,0)] + VFIELD[LM(i,j,k,0)]), Mesh[LM(i,j-2,k,1)], VFIELD[LM(i,j-2,k,0)], Mesh[LM(i,j-1,k,1)], VFIELD[LM(i,j-1,k,0)], Mesh[LM(i,j,k,1)], VFIELD[LM(i,j,k,0)], Mesh[LM(i,j+1,k,1)], VFIELD[LM(i,j+1,k,0)]);
            }
        }
    }

    ApplyBoundaries(WallProperty2, Bound2);  

    // MW Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 1; k < NZ; k++){
                WallProperty3.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j,k-1,2)] + Mesh[LM(i,j,k,2)]), 0.50 * (WFIELD[LM(i,j,k-1,0)] + WFIELD[LM(i,j,k,0)]), Mesh[LM(i,j,k-2,2)], WFIELD[LM(i,j,k-2,0)], Mesh[LM(i,j,k-2,2)], WFIELD[LM(i,j-1,k,0)], Mesh[LM(i,j,k,2)], WFIELD[LM(i,j,k,0)], Mesh[LM(i,j,k+1,2)], WFIELD[LM(i,j,k+1,0)]);
            }
        }
    }

    ApplyBoundaries(WallProperty3, Bound3);  

}

// Function to compute the Convective Term of a Property
void CFD_Solver::Get_ConvectiveTermScalar(Mesher MESH, Walls &WallProperty, double *ConvectiveProperty){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                ConvectiveProperty[LM(i,j,k,0)] = 1.0/MESH.Vol[LM(i,j,k,0)] * (
                                                + MESH.Surf[LM(i,j,k,0)] * U_W.Wall_U[LMU(i+1,j,k,0)] * WallProperty.Wall_U[LMU(i+1,j,k,0)]
                                                - MESH.Surf[LM(i,j,k,0)] * U_W.Wall_U[LMU(i,j,k,0)] * WallProperty.Wall_U[LMU(i,j,k,0)]
                                                + MESH.Surf[LM(i,j,k,1)] * V_W.Wall_V[LMV(i,j+1,k,0)] * WallProperty.Wall_V[LMV(i,j+1,k,0)]
                                                - MESH.Surf[LM(i,j,k,1)] * V_W.Wall_V[LMV(i,j,k,0)] * WallProperty.Wall_V[LMV(i,j,k,0)]
                                                + MESH.Surf[LM(i,j,k,2)] * W_W.Wall_W[LMW(i,j,k+1,0)] * WallProperty.Wall_W[LMW(i,j,k+1,0)]
                                                - MESH.Surf[LM(i,j,k,2)] * W_W.Wall_W[LMW(i,j,k,0)] * WallProperty.Wall_W[LMW(i,j,k,0)]
                                                );
            }
        }
    }

}