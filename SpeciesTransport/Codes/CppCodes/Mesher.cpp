//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MESH CREATION CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"

using namespace std;

#define LM(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo) * ((i) - Ix + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo) + ((Fx - Ix + 2*Halo) * (NY + 2*Halo) * (NZ + 2*Halo)) * (dim)

Mesher::Mesher(Memory M1, ReadData R1, Parallel P1){
	
    // Data necessary

        // Geometry data
        Xdominio = R1.GeometryData[1];
        Ydominio = R1.GeometryData[2];
        Zdominio = R1.GeometryData[3];

        // Meshing data
        NX = R1.ProblemNumericalData[1];
        NY = R1.ProblemNumericalData[2];
        NZ = R1.ProblemNumericalData[3];

        // Parallel computing data
        Rank = P1.Rank;
        Procesos = P1.Procesos;
        Ix = P1.Ix;
        Fx = P1.Fx;

        Halo = 2;


}

void Mesher::AllocateMemory(Memory M1){

    // Mesh 
    M_Nodes = M1.AllocateDouble(Fx - Ix + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 3);

}

void Mesher::Get_Mesh(){
int i, j, k;

    // Uniform Meshing
    if (Rank != 0 && Rank != Procesos -1){

        for (i = Ix - Halo; i < Fx + Halo; i++){
            for(j =  - Halo; j < NY + Halo; j++){
                for(k = - Halo; k < NZ + Halo; k++){
                    M_Nodes[M(i,j,k,0)] = (i/NX) * Xdominio; // X coordinate
                }
            }
        }

    }
    
}