//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MESH CREATION CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"

using namespace std;

#define LM(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo) + ((Fx[Rango] - Ix[Rango] + 2*Halo) * (NY + 2*Halo) * (NZ + 2*Halo)) * (dim)
#define GM(i,j,k,dim) (NY + 2*Halo) * (NZ + 2*Halo) * ((i) + Halo) + (NZ + 2*Halo) * ((j) + Halo) + ((k) + Halo) + (NY + 2*Halo) * (NZ + 2*Halo) * (NX + 2*Halo) * (dim)

Mesher::Mesher(ReadData R1, Parallel P1){
	
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

        Halo = 2;

}

void Mesher::AllocateMemory(Memory M1){

    if (Rango == 0){
        M_NodesGlobal = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 3);
        Test_MeshGlobal = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    }

    // Mesh 
    M_Nodes = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 3);

    // Scalar field
    Test_Mesh = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

void Mesher::Get_WorksplitInfo(Memory M1, Parallel P1){
int i;

    //Parallel computing
    Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

    for (i = 0; i < Procesos; i++){
        Ix[i] = P1.Ix[i];
        Fx[i] = P1.Fx[i];
    }

}

void Mesher::Get_Mesh(){
int i, j, k;

    // Uniform Meshing  
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for(j =  - Halo; j < NY + Halo; j++){
            for(k = - Halo; k < NZ + Halo; k++){
                M_Nodes[LM(i,j,k,0)] = 0.50 * ((2.0 * i + 1)/(NX + 1)) * Xdominio; // X coordinate
                M_Nodes[LM(i,j,k,1)] = 0.50 * ((2.0 * j + 1)/(NY + 1)) * Ydominio; // Y coordinate
                M_Nodes[LM(i,j,k,2)] = 0.50 * ((2.0 * k + 1)/(NZ + 1)) * Zdominio; // Z coordinate
            }
        }
    }

}

void Mesher::Get_ZeroGlobalMesh(){
int i, j, k;

    // Uniform Meshing  
    for (i = - Halo; i < NX + Halo; i++){
        for(j =  - Halo; j < NY + Halo; j++){
            for(k = - Halo; k < NZ + Halo; k++){   
                M_NodesGlobal[GM(i,j,k,0)] = 0.50 * ((2.0 * i + 1)/(NX + 1)) * Xdominio; // X coordinate
                M_NodesGlobal[GM(i,j,k,1)] = 0.50 * ((2.0 * j + 1)/(NY + 1)) * Ydominio; // Y coordinate
                M_NodesGlobal[GM(i,j,k,2)] = 0.50 * ((2.0 * k + 1)/(NZ + 1)) * Zdominio; // Z coordinate
            }
        }
    }

}

void Mesher::RunMesher(Memory M1, Parallel P1){

    Get_WorksplitInfo(M1, P1);
    AllocateMemory(M1); 
    Get_Mesh();

    if (Rango == 0){
        Get_ZeroGlobalMesh();
    }

}