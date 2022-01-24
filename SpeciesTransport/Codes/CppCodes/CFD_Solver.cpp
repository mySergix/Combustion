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

CFD_Solver::CFD_Solver(Memory M1, ReadData R1, Parallel P1, Species_Solver SPE_S1){

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
//#include "CFD_Solver_Utilities.cpp"
//#include "CFD_Solver_BoundaryConditions.cpp"
//#include "CFD_Solver_Mass.cpp"
//#include "CFD_Solver_Momentum.cpp"
//#include "CFD_Solver_Energy.cpp"

// Function to communicate all the velocity fields (local fields)
void CFD_Solver::CommunicateVelocities(Parallel P1, double *UFIELD, double *VFIELD, double *WFIELD){

    P1.CommunicateLocalMatrix(UFIELD, UFIELD); // Velocity U
    P1.CommunicateLocalMatrix(VFIELD, VFIELD); // Velocity V
    P1.CommunicateLocalMatrix(WFIELD, WFIELD); // Velocity W

}