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

        struct Global
        {
            double *Density;
            double *U;
            double *V;
            double *W;
        };

		//Constructor de la clase
		PostProcess(Memory, ReadData, Parallel);
		
		//Metodos de la clase
        void GlobalEscalarVTK(Mesher, string, string, string, double*, int);
        void GlobalVectorialVTK(Mesher, string, string, string, Global&, int);
        
};