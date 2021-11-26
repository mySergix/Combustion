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

Mesher::Mesher(Memory M1, ReadData R1, Parallel P1){
	
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

}

// Function to allocate the memory of the mesher
void Mesher::AllocateMemory(Memory M1){

    if (Rango == 0){
        // Global mesh (core 0)
        Global_Node_Mesh = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 3);
    }

    // Local meshes
    Node_Mesh = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 3);
    DeltaP = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 3); // Deltas (distance) of the nodal mesh
    Surf = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 3); // Finite volume walls surface
    Vol = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1); // Finite volumes volume

}

// Function to create the nodal local meshes
void Mesher::Get_LocalMeshes(){
int i, j, k;

    // Uniform Meshing  
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for(j =  - Halo; j < NY + Halo; j++){
            for(k = - Halo; k < NZ + Halo; k++){
                Node_Mesh[LM(i,j,k,0)] = 0.50 * ((2.0 * i + 1)/(NX + 1)) * Xdominio; // X coordinate
                Node_Mesh[LM(i,j,k,1)] = 0.50 * ((2.0 * j + 1)/(NY + 1)) * Ydominio; // Y coordinate
                Node_Mesh[LM(i,j,k,2)] = 0.50 * ((2.0 * k + 1)/(NZ + 1)) * Zdominio; // Z coordinate
            }
        }
    }

}

// Function to calculate the deltas of the local meshes
void Mesher::Get_LocalMesh_Deltas(){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                DeltaP[LM(i,j,k,0)] = 0.50 * (Node_Mesh[LM(i+1,j,k,0)] - Node_Mesh[LM(i-1,j,k,0)]);
                DeltaP[LM(i,j,k,1)] = 0.50 * (Node_Mesh[LM(i,j+1,k,1)] - Node_Mesh[LM(i,j-1,k,1)]);
                DeltaP[LM(i,j,k,2)] = 0.50 * (Node_Mesh[LM(i,j,k+1,2)] - Node_Mesh[LM(i,j,k-1,2)]);
            }
        }
    }
}

// Function to calculate the surfaces of the control volumes of the mesh
void Mesher::Get_LocalMesh_Surfaces(){
int i, j, k;
 
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Surf[LM(i,j,k,0)] = DeltaP[LM(i,j,k,1)] * DeltaP[LM(i,j,k,2)];
                Surf[LM(i,j,k,1)] = DeltaP[LM(i,j,k,0)] * DeltaP[LM(i,j,k,2)];
                Surf[LM(i,j,k,2)] = DeltaP[LM(i,j,k,0)] * DeltaP[LM(i,j,k,1)];
            }
        }
    }

}

// Function to calculate the volume of the control volumes of the mesh
void Mesher::Get_LocalMesh_Volumes(){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Vol[LM(i,j,k,0)] = DeltaP[LM(i,j,k,0)] * DeltaP[LM(i,j,k,1)] * DeltaP[LM(i,j,k,2)];
            }
        }
    }

}

// Function to calculate the global mesh for core 0
void Mesher::Get_GlobalMesh(){
int i, j, k;

    // Uniform Meshing  
    for (i = - Halo; i < NX + Halo; i++){
        for(j =  - Halo; j < NY + Halo; j++){
            for(k = - Halo; k < NZ + Halo; k++){   
                Global_Node_Mesh[GM(i,j,k,0)] = 0.50 * ((2.0 * i + 1)/(NX + 1)) * Xdominio; // X coordinate
                Global_Node_Mesh[GM(i,j,k,1)] = 0.50 * ((2.0 * j + 1)/(NY + 1)) * Ydominio; // Y coordinate
                Global_Node_Mesh[GM(i,j,k,2)] = 0.50 * ((2.0 * k + 1)/(NZ + 1)) * Zdominio; // Z coordinate
            }
        }
    }

}

// Function to run all the mesher and get all the necessary matrix
void Mesher::RunMesher(Memory M1){

    AllocateMemory(M1); 
    Get_LocalMeshes();
    Get_LocalMesh_Deltas();
    Get_LocalMesh_Surfaces();
    Get_LocalMesh_Volumes();

    if (Rango == 0){
        Get_GlobalMesh();
    }

    
}