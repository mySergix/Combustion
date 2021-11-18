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
        int Rank;
        int Procesos;

        int Halo;

		//Constructor de la clase
		Mesher(Memory M1, ReadData R1, Parallel P1);
		
		//Metodos de la clase
        void AllocateMemory(Memory);
	
			
};