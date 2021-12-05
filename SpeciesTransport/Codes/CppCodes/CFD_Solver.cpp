//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS                                      //
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

        Beta = 0.6;
        GlobalConvergence = 1e-5;

        c = 340;

};

// Function to allocate memory for all matrix
void CFD_Solver::AllocateMemory(Memory M1){

    // Density
    Allocate_StructureMemory(M1, Density);

    // Velocity U
    Allocate_StructureMemory(M1, U);

    // Velocity V
    Allocate_StructureMemory(M1, V);

    // Velocity W
    Allocate_StructureMemory(M1, W);

    if (Rango == 0){
        GlobalMatrix.Density = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

        GlobalMatrix.U = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
        GlobalMatrix.V = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
        GlobalMatrix.W = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    }

    Pressure.Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    Pressure.Gradient_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Pressure.Gradient_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Pressure.Gradient_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    // Viscous Stresses
    Stress.Divergence = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_xx = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1); 
    Stress.Tau_yy = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_zz = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    Stress.Tau_xy_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_xy_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);

    Stress.Tau_yz_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_yz_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    Stress.Tau_zx_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);
    Stress.Tau_zx_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

// Function to allocate memory for a certain structure
void CFD_Solver::Allocate_StructureMemory(Memory M1, Property &PropertyName){

    PropertyName.Past = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Fut = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    PropertyName.Wall_U = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1); 
    PropertyName.Wall_V = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Wall_W = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    PropertyName.Diff_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Diff_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Diff_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    PropertyName.Convective = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    PropertyName.ContributionPast = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.ContributionPres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    
    PropertyName.Bottom = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, 1, NZ + 2*Halo, 1);
    PropertyName.Top = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, 1, NZ + 2*Halo, 1);

    PropertyName.Here = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, 1, 1);
    PropertyName.There = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, 1, 1);

    if (Rango == 0){
        PropertyName.Left = M1.AllocateDouble(1, NY + 2*Halo, NZ + 2*Halo, 1);     
    }
    if (Rango == Procesos - 1){
        PropertyName.Right = M1.AllocateDouble(1, NY + 2*Halo, NZ + 2*Halo, 1);  
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
    PhiC_ = (PHIC - PHIU) / (PHID - PHIU + 1e-10);
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
            Density.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // Density
            U.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // U Velocity
            V.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // V Velocity
            W.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // V Velocity

            //Parte Top
            Density.Top[TOP_BOT(i,NY,k,0)] = 0.0; //Density
            U.Top[TOP_BOT(i,NY,k,0)] = 0.0; // U Velocity
            V.Top[TOP_BOT(i,NY,k,0)] = 0.0; // V Velocity
            W.Top[TOP_BOT(i,NY,k,0)] = 0.0; // W Velocity

        }
    }

    // Parte Here y There
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){

            // Parte Here
            Density.Here[HER_THE(i,j,0,0)] = 0.0; // Density
            

            // Parte There 
            Density.There[HER_THE(i,j,NZ,0)] = 0.0; // Density
            U.There[HER_THE(i,j,NZ,0)] = 0.0; // Velocity U
            V.There[HER_THE(i,j,NZ,0)] = 0.0; // Velocity V
            W.There[HER_THE(i,j,NZ,0)] = 0.0; // Velocity W

        }
    }

    // Parte Left
    if (Rango == 0){

        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){

                Density.Left[LEF_RIG(0,j,k,0)] = 1.0; // Density
                U.Left[LEF_RIG(0,j,k,0)] = 10.0; // Velocity U
                V.Left[LEF_RIG(0,j,k,0)] = 0.0; // Velocity V
                W.Left[LEF_RIG(0,j,k,0)] = 0.0; // Velocity W

            }
        }

    }

}

// Function to calculate and update periodic conditions
void CFD_Solver::Get_PeriodicConditions(Mesher MESH, Property &PropertyName){
int i, j, k;

    // Here
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < 0; k++){
                PropertyName.Pres[LM(i,j,k,0)] = PropertyName.Pres[LM(i,j,NZ+k,0)];
            }
        }
    }

    // There
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = NZ; k < NZ + Halo; k++){
                PropertyName.Pres[LM(i,j,k,0)] = PropertyName.Pres[LM(i,j,k-NZ,0)];
            }
        }
    } 

    // Walls
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            PropertyName.Here[HER_THE(i,j,0,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,-1,2)] + MESH.Node_Mesh[LM(i,j,0,2)]), 0.50 * (W.Pres[LM(i,j,-1,0)] + W.Pres[LM(i,j,0,0)]), MESH.Node_Mesh[LM(i,j,-2,2)], PropertyName.Pres[LM(i,j,-2,0)], MESH.Node_Mesh[LM(i,j,-1,2)], PropertyName.Pres[LM(i,j,-1,0)], MESH.Node_Mesh[LM(i,j,0,2)], PropertyName.Pres[LM(i,j,0,0)], MESH.Node_Mesh[LM(i,j,1,2)], PropertyName.Pres[LM(i,j,1,0)]);
            PropertyName.There[HER_THE(i,j,NZ,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,NZ-1,2)] + MESH.Node_Mesh[LM(i,j,NZ,2)]), 0.50 * (W.Pres[LM(i,j,NZ-1,0)] + W.Pres[LM(i,j,NZ,0)]), MESH.Node_Mesh[LM(i,j,NZ-2,2)], PropertyName.Pres[LM(i,j,NZ-2,0)], MESH.Node_Mesh[LM(i,j,NZ-1,2)], PropertyName.Pres[LM(i,j,NZ-1,0)], MESH.Node_Mesh[LM(i,j,NZ,2)], PropertyName.Pres[LM(i,j,NZ,0)], MESH.Node_Mesh[LM(i,j,NZ+1,2)], PropertyName.Pres[LM(i,j,NZ+1,0)]);
        }
    }

}

// Function to update continuously the boundary conditions
void CFD_Solver::Update_BoundaryConditions(Mesher MESH){
int i, j, k;

    Get_PeriodicConditions(MESH, U); // Periodic Conditions Velocity U
    Get_PeriodicConditions(MESH, V); // Periodic Conditions Velocity V
    Get_PeriodicConditions(MESH, W); // Periodic Conditions Velocity W

    if (Rango == Procesos - 1){

        // Right Side
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Density.Right[LEF_RIG(0,j,k,0)] = Density.Pres[LM(NX-1,j,k,0)]; // Density
                U.Right[LEF_RIG(0,j,k,0)] = U.Pres[LM(NX-1,j,k,0)]; // Velocity U
                V.Right[LEF_RIG(0,j,k,0)] = V.Pres[LM(NX-1,j,k,0)]; // Velocity V
                W.Right[LEF_RIG(0,j,k,0)] = W.Pres[LM(NX-1,j,k,0)]; // Velocity W
            }
        }
        
        // Pressure Right Side
        for (i = NX; i < NX + Halo; i++){
            for (j = - Halo; j < NY + Halo; j++){
                for (k = - Halo; k < NZ + Halo; k++){
                    Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(NX-1,j,k,0)]; // Pressure
                }
            }
        }

    }

    if (Rango == 0){

        // Pressure Left Side
        for (i = - Halo; i < 0; i++){
            for (j = - Halo; j < NY + Halo; j++){
                for (k = - Halo; k < NZ + Halo; k++){
                    Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(0,j,k,0)]; // Pressure
                }
            }
        }

    }

    // Bottom 
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < 0; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(i,0,k,0)];
            }
        }
    }

    // Pressure Top
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(i,NY-1,k,0)];
            }
        }
    }

    // Pressure Here
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < 0; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[i,j,NX + k,0];
            }
        }
    }

    // Pressure There
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = NX; k < NX + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[i,j,k - NX,0];
            }
        }
    }

}

// Function to set initial values for the variables and properties fields
void CFD_Solver::Set_InitialValues(){
int i, j, k;

    // Densities
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Density.Pres[LM(i,j,k,0)] = 0.95;
                Density.Past[LM(i,j,k,0)] = 0.95;
            }
        }
    }

    // Velocities
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){   
            for (k = 0; k < NZ; k++){

                // Velocity U
                U.Pres[LM(i,j,k,0)] = 10.0;
                U.Past[LM(i,j,k,0)] = 10.0;

                // Velocity V
                V.Pres[LM(i,j,k,0)] = 0.0;
                V.Past[LM(i,j,k,0)] = 0.0;

                // Velocity W
                W.Pres[LM(i,j,k,0)] = 0.0;
                W.Past[LM(i,j,k,0)] = 0.0;
            }
        }
    }

    // Pressure
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] =  Density.Pres[LM(i,j,k,0)] * 287.0 * 298.15;
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

// Function to calulate the time step of the simulation
void CFD_Solver::Get_TimeStep(Mesher MESH){
int i, j, k;
double Courant = 0.2;

DeltaT = 1000.0;

MPI_Status ST;

    // Comparacion
    // Delta += (NewDelta - Delta) * (NewDelta < Delta)

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){

                // CFL Convective

                    // Velocity U
                    DeltaT += (Courant * (MESH.DeltaP[LM(i,j,k,0)] / (abs(U.Pres[LM(i,j,k,0)]) + c + 1e-10)) - DeltaT) * (Courant * (MESH.DeltaP[LM(i,j,k,0)] / (abs(U.Pres[LM(i,j,k,0)]) + c + 1e-10)) < DeltaT);

                    // Velocity V
                    DeltaT += (Courant * (MESH.DeltaP[LM(i,j,k,1)] / (abs(V.Pres[LM(i,j,k,0)]) + c + 1e-10)) - DeltaT) * (Courant * (MESH.DeltaP[LM(i,j,k,1)] / (abs(V.Pres[LM(i,j,k,0)]) + c + 1e-10)) < DeltaT);

                    // Velocity W
                    DeltaT += (Courant * (MESH.DeltaP[LM(i,j,k,2)] / (abs(W.Pres[LM(i,j,k,0)]) + c + 1e-10)) - DeltaT) * (Courant * (MESH.DeltaP[LM(i,j,k,2)] / (abs(W.Pres[LM(i,j,k,0)]) + c + 1e-10)) < DeltaT);
            
            }
        }
    }

    MPI_Allreduce(&DeltaT, &DeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

// Calculation of Property Value at the Walls of the Volumes
void CFD_Solver::Get_WallsValue_Scalar(Parallel P1, double *Mesh, Property &PropertyName){
int i, j, k;

    // Communication of the local metrix of the properties
    P1.CommunicateLocalMatrix(PropertyName.Pres, PropertyName.Pres);

    // MU Walls
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (Mesh[LM(i - 1,j,k,0)] + Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), Mesh[LM(i-2,j,k,0)], PropertyName.Pres[LM(i-2,j,k,0)], Mesh[LM(i-1,j,k,0)], PropertyName.Pres[LM(i-1,j,k,0)], Mesh[LM(i,j,k,0)], PropertyName.Pres[LM(i,j,k,0)], Mesh[LM(i+1,j,k,0)], PropertyName.Pres[LM(i+1,j,k,0)]); 
            }
        }
    }

    // MV Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j-1,k,1)] + Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), Mesh[LM(i,j-2,k,1)], PropertyName.Pres[LM(i,j-2,k,0)], Mesh[LM(i,j-1,k,1)], PropertyName.Pres[LM(i,j-1,k,0)], Mesh[LM(i,j,k,1)], PropertyName.Pres[LM(i,j,k,0)], Mesh[LM(i,j+1,k,1)], PropertyName.Pres[LM(i,j+1,k,0)]);       
            }
        }
    }

    // MW Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 1; k < NZ; k++){
                PropertyName.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j,k-1,2)] + Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), Mesh[LM(i,j,k-2,2)], PropertyName.Pres[LM(i,j,k-2,0)], Mesh[LM(i,j,k-1,2)], PropertyName.Pres[LM(i,j,k-1,0)], Mesh[LM(i,j,k,2)], PropertyName.Pres[LM(i,j,k,0)], Mesh[LM(i,j,k+1,2)], PropertyName.Pres[LM(i,j,k+1,0)]);
            }
        }
    }

    ApplyBoundaries(PropertyName);   

}

// Function for the application of the boundary condition to the walls of the volumes
void CFD_Solver::ApplyBoundaries(Property &PropertyName){
int i, j, k;

    // Left Side
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){               
                PropertyName.Wall_U[LMU(0,j,k,0)] = PropertyName.Left[LEF_RIG(0,j,k,0)];
            }
        }
    }

    // Right Side
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){               
                PropertyName.Wall_U[LMU(NX,j,k,0)] = PropertyName.Right[LEF_RIG(NX,j,k,0)];
            }
        }
    }

    // Top and Bottom Sides
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            // Bottom Part
            PropertyName.Wall_V[LMV(i,0,k,0)] = PropertyName.Bottom[TOP_BOT(i,0,k,0)];

            // Top Part
            PropertyName.Wall_V[LMV(i,NY,k,0)] = PropertyName.Top[TOP_BOT(i,NY,k,0)];
        }
    }

    // Here and There Sides
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            // Here Part
            PropertyName.Wall_W[LMW(i,j,0,0)] = PropertyName.Here[HER_THE(i,j,0,0)];

            // There Part
            PropertyName.Wall_W[LMW(i,j,NZ,0)] = PropertyName.There[HER_THE(i,j,NZ,0)];
        }
    }

}

// Calculation of Velocities at the Walls of the Volumes
void CFD_Solver::Get_WallsVelocities(Mesher MESH){
int i, j, k;

    // MU Walls
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], U.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], U.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], U.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], U.Pres[LM(i+1,j,k,0)]);
                V.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], V.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], V.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], V.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], V.Pres[LM(i+1,j,k,0)]);
                W.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], W.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], W.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], W.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], W.Pres[LM(i+1,j,k,0)]);
            }
        }
    }

    ApplyBoundaries(U);  

    // MV Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], U.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], U.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], U.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], U.Pres[LM(i,j+1,k,0)]);
                V.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], V.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], V.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], V.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], V.Pres[LM(i,j+1,k,0)]);
                W.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], W.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], W.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], W.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], W.Pres[LM(i,j+1,k,0)]);
            }
        }
    }

    ApplyBoundaries(V);  

    // MW Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 1; k < NZ; k++){
                U.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], U.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], U.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], U.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], U.Pres[LM(i,j,k+1,0)]);
                V.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], V.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], V.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], V.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], V.Pres[LM(i,j,k+1,0)]);
                W.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], W.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], W.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], W.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], W.Pres[LM(i,j,k+1,0)]);
            }
        }
    }

    ApplyBoundaries(W);  

}

// Function to compute the Convective Term of a Property
void CFD_Solver::Get_ConvectiveTerm(Mesher MESH, Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Convective[LM(i,j,k,0)] = 1.0/MESH.Vol[LM(i,j,k,0)] * (
                                                     + MESH.Surf[LM(i,j,k,0)] * U.Wall_U[LMU(i+1,j,k,0)] * PropertyName.Wall_U[LMU(i+1,j,k,0)]
                                                     - MESH.Surf[LM(i,j,k,0)] * U.Wall_U[LMU(i,j,k,0)] * PropertyName.Wall_U[LMU(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,1)] * V.Wall_V[LMV(i,j+1,k,0)] * PropertyName.Wall_V[LMV(i,j+1,k,0)]
                                                     - MESH.Surf[LM(i,j,k,1)] * V.Wall_V[LMV(i,j,k,0)] * PropertyName.Wall_V[LMV(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,2)] * W.Wall_W[LMW(i,j,k+1,0)] * PropertyName.Wall_W[LMW(i,j,k+1,0)]
                                                     - MESH.Surf[LM(i,j,k,2)] * W.Wall_W[LMW(i,j,k,0)] * PropertyName.Wall_W[LMW(i,j,k,0)]
                                                     );
            }
        }
    }

}

//Function to calculate the pressure
void CFD_Solver::Get_Pressure(){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Pressure.Pres[LM(i,j,k,0)] = Density.Pres[LM(i,j,k,0)] * 287.0 * 298.15;
            }
        }
    }

}

// Function to calculate the pressure gradient
void CFD_Solver::Get_PressureGradient(Mesher MESH){
int i, j, k;
double Px_1, Px_2, Py_1, Py_2, Pz_1, Pz_2;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Px_1 = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], Pressure.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], Pressure.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], Pressure.Pres[LM(i+1,j,k,0)]);
                Px_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,0)] + MESH.Node_Mesh[LM(i + 1,j,k,0)]), 0.50 * (U.Pres[LM(i,j,k,0)] + U.Pres[LM(i + 1,j,k,0)]), MESH.Node_Mesh[LM(i-1,j,k,0)], Pressure.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], Pressure.Pres[LM(i+1,j,k,0)], MESH.Node_Mesh[LM(i+2,j,k,0)], Pressure.Pres[LM(i+2,j,k,0)]);
                
                Py_1 = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], Pressure.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], Pressure.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], Pressure.Pres[LM(i,j+1,k,0)]);
                Py_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,1)] + MESH.Node_Mesh[LM(i,j+1,k,1)]), 0.50 * (V.Pres[LM(i,j,k,0)] + V.Pres[LM(i,j+1,k,0)]), MESH.Node_Mesh[LM(i,j-1,k,1)], Pressure.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], Pressure.Pres[LM(i,j+1,k,0)], MESH.Node_Mesh[LM(i,j+2,k,1)], Pressure.Pres[LM(i,j+2,k,0)]);
                
                Pz_1 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], Pressure.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], Pressure.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], Pressure.Pres[LM(i,j,k+1,0)]);
                Pz_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,2)] + MESH.Node_Mesh[LM(i,j,k+1,2)]), 0.50 * (W.Pres[LM(i,j,k,0)] + W.Pres[LM(i,j,k+1,0)]), MESH.Node_Mesh[LM(i,j,k-1,2)], Pressure.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], Pressure.Pres[LM(i,j,k+1,0)], MESH.Node_Mesh[LM(i,j,k+2,2)], Pressure.Pres[LM(i,j,k+2,0)]);

                Pressure.Gradient_X[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,0)]) * (Px_2 - Px_1);
                Pressure.Gradient_Y[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,1)]) * (Py_2 - Py_1);
                Pressure.Gradient_Z[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,2)]) * (Pz_2 - Pz_1);
            }
        }
    }

}

// Function to calculate the divergence of the field
void CFD_Solver::Get_Divergence(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Divergence[LM(i,j,k,0)] = (1.0 / VolMP[LM(i,j,k,0)]) * (
                                        - MESH.Surf[LM(i,j,k,0)] * U.Wall_U[LMU(i,j,k,0)]
                                        + MESH.Surf[LM(i,j,k,0)] * U.Wall_U[LMU(i+1,j,k,0)]
                                        - MESH.Surf[LM(i,j,k,1)] * V.Wall_V[LMV(i,j,k,0)]
                                        + MESH.Surf[LM(i,j,k,1)] * V.Wall_V[LMV(i,j+1,k,0)]
                                        - MESH.Surf[LM(i,j,k,2)] * W.Wall_W[LMW(i,j,k,0)]
                                        + MESH.Surf[LM(i,j,k,2)] * W.Wall_W[LMW(i,j,k+1,0)]
                                        );
            }
        }
    }

}

// Function to calculate the velociy gradients at the walls
void CFD_Solver::Get_VelocityGradients(Mesher MESH, Property &PropertyName){
int i, j, k;

    // Gradients in X direction
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Diff_X[LMU(i,j,k,0)] = (PropertyName.Pres[LM(i,j,k,0)] - PropertyName.Pres[LM(i-1,j,k,0)]) / (0.50 * (MESH.DeltaP[LM(i,j,k,0)] + MESH.DeltaP[LM(i-1,j,k,0)]));
            }
        }
    }

    // Gradients in Y direction
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                // Parte Top
                PropertyName.Diff_Y[LMV(i,0,k,0)] = 0.0;

                // Parte Bottom
                PropertyName.Diff_Y[LMV(i,NY,k,0)] = 0.0;

                // Resto
                PropertyName.Diff_Y[LMV(i,j,k,0)] = (PropertyName.Pres[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j-1,k,0)]) / (0.50 * (MESH.DeltaP[LM(i,j,k,1)] + MESH.DeltaP[LM(i,j-1,k,1)]));
            }
        }
    }

    // Gradients in Z direction
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                PropertyName.Diff_Z[LMW(i,j,k,0)] = (PropertyName.Pres[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j,k-1,0)]) / (0.50 * (MESH.DeltaP[LM(i,j,k,2)] + MESH.DeltaP[LM(i,j,k-1,2)]));
            }
        }
    }

}

// Function to calculate the viscous stresses of the fluid
void CFD_Solver::Get_ViscousStresses(Mesher MESH){
int i, j, k;

    // Tau XX
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_xx[LMU(i,j,k,0)] =  - (2.0/3.0) * mu * 0.50 * (Divergence[LM(i-1,j,k,0)] + Divergence[LM(i,j,k,0)]) + 2.0 * mu * U.Diff_X[LMU(i,j,k,0)];
            }
        }
    }

    // Tau YY
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY + 1; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_yy[LMV(i,j,k,0)] =  - (2.0/3.0) * mu * 0.50 * (Divergence[LM(i,j-1,k,0)] + Divergence[LM(i,j,k,0)]) + 2.0 * mu * V.Diff_Y[LMV(i,j,k,0)];
            }
        }
    }

    // Tau ZZ
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                Stress.Tau_zz[LMW(i,j,k,0)] =  - (2.0/3.0) * mu * 0.50 * (Divergence[LM(i,j,k-1,0)] + Divergence[LM(i,j,k,0)]) + 2.0 * mu * W.Diff_Z[LMW(i,j,k,0)];
            }
        }
    }

    // Tau XY
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_xy_Y[LMU(i,j,k,0)] =  mu * (U.Diff_Y[LMV(i,j,k,0)] + ;
            }
        }
    }

}
// Function to calculate the step contribution to the density
void CFD_Solver::Get_StepContribution_Density(Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.ContributionPres[LM(i,j,k,0)] = - PropertyName.Convective[LM(i,j,k,0)];
            }
        }
    }

}

// Function to calculate the step contribution to each velocity
void CFD_Solver::Get_StepContribution_Velocity(Property &PropertyName, double *PressureGradient){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.ContributionPres[LM(i,j,k,0)] = (1.0 / Density.Pres[LM(i,j,k,0)]) * (- PropertyName.Convective[LM(i,j,k,0)] - PressureGradient[LM(i,j,k,0)]) + PropertyName.Gravity;
            }
        }
    }

}

// Function to integrate the density calculation
void CFD_Solver::Get_TemporalIntegration(Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Fut[LM(i,j,k,0)] = (2.0 * Beta * PropertyName.Pres[LM(i,j,k,0)] - (Beta - 0.50) * PropertyName.Past[LM(i,j,k,0)]) / (Beta + 0.50) 
                                              + DeltaT * ((1.0 + Beta) * PropertyName.ContributionPres[LM(i,j,k,0)] - Beta * PropertyName.ContributionPast[LM(i,j,k,0)]); 
            }
        }
    }

}

// Function to update the values of the properties fields
void CFD_Solver::UpdateField(Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Past[LM(i,j,k,0)] = PropertyName.Pres[LM(i,j,k,0)];
                PropertyName.Pres[LM(i,j,k,0)] = PropertyName.Fut[LM(i,j,k,0)];

                PropertyName.ContributionPast[LM(i,j,k,0)] = PropertyName.ContributionPres[LM(i,j,k,0)];
            }
        }
    }

}

// Function to calculate the difference between time steps
void CFD_Solver::Get_MaximumDifference(Property &PropertyName, double &MaxDiff){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                MaxDiff += (abs((PropertyName.Fut[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j,k,0)]) / (PropertyName.Pres[LM(i,j,k,0)] + 1e-10)) - MaxDiff) * (abs((PropertyName.Fut[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j,k,0)]) / (PropertyName.Pres[LM(i,j,k,0)] + 1e-10)) > MaxDiff);
            }
        }
    }

    MPI_Allreduce(&MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

}

// Function to Check for convergence criteria
void CFD_Solver::Get_ConvergenceCriteria(){
int i, j, k;
MaxDiff = 0.0;

    // Check Convergence for Density
    Get_MaximumDifference(Density, MaxDiff);

}

// Function to run the CFD Solver
void CFD_Solver::RunSolver(Memory M1, Parallel P1, Mesher MESH, PostProcess PP1){

    double Time = 0.0; // Time of the simulation
    int Step = 0; // Steps of the simulation

    MaxDiff = 2.0 * GlobalConvergence;

    char FileName_1[300];

    AllocateMemory(M1); // Allocate memory for the matrix
    Set_InitialValues(); // Get the inital values of the properties
    Get_BoundaryConditions(); // Get the fixed boundary conditions
    
    
    while (MaxDiff >= GlobalConvergence){

        Step++;
        Update_BoundaryConditions(MESH);
        Get_TimeStep(MESH);
        Time += DeltaT;

        CommunicateVelocities(P1, U.Pres, V.Pres, W.Pres);
        Get_WallsVelocities(MESH);
        Get_WallsValue_Scalar(P1, MESH.Node_Mesh, Density);

        // Density
        Get_ConvectiveTerm(MESH, Density);
        Get_StepContribution_Density(Density);
        Get_TemporalIntegration(Density);

        // Velocities
        //Get_Pressure();
        //Get_PressureGradient(MESH);

        // Velocity U
        Get_ConvectiveTerm(MESH, U);
        Get_StepContribution_Velocity(U, Pressure.Gradient_X);

        // Velocity V
        Get_ConvectiveTerm(MESH, V);
        Get_StepContribution_Velocity(V, Pressure.Gradient_Y);

        // Velocity W
        Get_ConvectiveTerm(MESH, W);
        Get_StepContribution_Velocity(W, Pressure.Gradient_Z);

        Get_TemporalIntegration(U);
        Get_TemporalIntegration(V);
        Get_TemporalIntegration(W);

        if (Step % 100 == 0){
            Get_ConvergenceCriteria();
            P1.SendMatrixToZero(Density.Pres, GlobalMatrix.Density);
            
            P1.SendMatrixToZero(U.Pres, GlobalMatrix.U);
            P1.SendMatrixToZero(V.Pres, GlobalMatrix.V);
            P1.SendMatrixToZero(W.Pres, GlobalMatrix.W);

            if (Rango == 0){
                cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiff<<endl;
                sprintf(FileName_1, "Density_Field_Step_%d", Step);
                PP1.GlobalEscalarVTK(MESH, "CombustionResults/", "Densidad", FileName_1, GlobalMatrix.Density, 0);

                sprintf(FileName_1, "Velocities_Field_Step_%d", Step);
                PP1.GlobalVectorialVTK(MESH, "CombustionResults/", "Velocities", FileName_1, GlobalMatrix.U, GlobalMatrix.V, GlobalMatrix.W, 0);
            }

            MPI_Barrier(MPI_COMM_WORLD);	
        }
        
        UpdateField(Density);
        UpdateField(U);
        UpdateField(V);
        UpdateField(W);

    }
    
}
