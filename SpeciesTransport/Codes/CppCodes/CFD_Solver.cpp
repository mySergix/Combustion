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
        mu = 1.0;

};

// Parts of the CFD SOLVER Class
#include "CFD_Solver_Memory.cpp"
#include "CFD_Solver_Utilities.cpp"
#include "CFD_Solver_BoundaryConditions.cpp"
#include "CFD_Solver_Mass.cpp"
#include "CFD_Solver_Momentum.cpp"
#include "CFD_Solver_Energy.cpp"

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

  
}



// Function to communicate all the velocity fields (local fields)
void CFD_Solver::CommunicateVelocities(Parallel P1, double *UFIELD, double *VFIELD, double *WFIELD){

    P1.CommunicateLocalMatrix(UFIELD, UFIELD); // Velocity U
    P1.CommunicateLocalMatrix(VFIELD, VFIELD); // Velocity V
    P1.CommunicateLocalMatrix(WFIELD, WFIELD); // Velocity W

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
