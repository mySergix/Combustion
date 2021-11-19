//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR POST PROCESSING CLASS                              //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

class PostProcess{	

	public:

        // Variables de la clase

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

		//Constructor de la clase
		PostProcess(Memory M1, ReadData R1, Parallel P1);
		
		//Metodos de la clase
        void GlobalEscalarVTK(string, string, string, double*, double*, int);
        	
};