//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR PARALLEL PROGRAMMING CLASS                         //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include </usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h>

using namespace std;

class Parallel{	

	public:

        // Variables de la clase

            // Meshing data
            int NX;
            int NY;
            int NZ;

            // Parallel computing data
		    int Halo;
            int Rango;
            int Procesos;

            int* Ix;
            int* Fx;
            


		//Constructor de la clase
		Parallel(Memory M1);
		
		//Metodos de la clase
        void Rango_Procesos();
        void Total_Procesos();
        void AllocateMemory(Memory);
        void WorkSplit(int, int&, int&);
			
};