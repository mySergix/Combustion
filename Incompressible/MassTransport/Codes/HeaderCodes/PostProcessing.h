#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

class PostProcessing{	

	private:

	public:
		//Constructor de la clase
		PostProcessing(Memory, ReadData, Mesher);
		
		//Datos de la clase
		int Problema;

    	int NX;
    	int NY;
    	int NZ;

		int Halo;
		int HP;
		
		//Metodos de la clase
		void VTK_GlobalScalar3D(string, string, string, Mesher, double*);
		void VTK_GlobalVectorial3D(string, string, string, Mesher, double*, double*, double*);
		
};
