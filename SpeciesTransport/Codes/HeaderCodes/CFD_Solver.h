//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR CFD SOLVER CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

class CFD_Solver{	

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
        double* Test_Mesh;
        double* Test_MeshGlobal;

		//Constructor de la clase
		CFD_Solver(Memory M1, ReadData R1, Parallel P1);
		
		//Metodos de la clase
        
			
};