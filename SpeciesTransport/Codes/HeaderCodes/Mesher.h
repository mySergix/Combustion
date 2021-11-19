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
        double* M_Nodes;
        double* M_NodesGlobal;

        double* Test_Mesh;
        double* Test_MeshGlobal;

		//Constructor de la clase
		Mesher(ReadData R1, Parallel P1);
		
		//Metodos de la clase
        void AllocateMemory(Memory);
        void Get_WorksplitInfo(Memory, Parallel);
        void Get_Mesh();
        void Get_ZeroGlobalMesh();
        void RunMesher(Memory, Parallel);
			
};