//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR MESH CREATION CLASS                                //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

class Mesher{	

	public:

        // Variables de la clase

        // Geometry data
        double Xdominio;
        double Ydominio;
        double Zdominio;

        // Meshing data
        int NX;
        int NY;
        int NZ; 

        // Parallel computing data
        int Rango;
        int Procesos;
        int* Ix;
        int* Fx;

        int Halo;

        //Matrices necesarias
        double *Node_Mesh; // Nodal mesh
        double *DeltaP; // Deltas (distance) of the nodal mesh
        double *Surf; // Finite volume walls surface
        double *Vol; // Finite volumes volume
        
        // Global matrix
        double *Global_Node_Mesh;

		//Constructor de la clase
		Mesher(Memory, ReadData, Parallel);
		
		//Metodos de la clase
        void AllocateMemory(Memory);
        
        void Get_LocalMeshes();
        void Get_LocalMesh_Deltas();
        void Get_LocalMesh_Surfaces();
        void Get_LocalMesh_Volumes();

        void Get_GlobalMesh();

        void RunMesher(Memory);
			
};