//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR PARALLEL PROGRAMMING CLASS                         //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

class Parallel{
	private:
		
		
	public:

		//Datos generales del problema
		int NX;
		int NY;
		int NZ;
		
		//Datos y variables de computación paralela
		int Rango;
		int Procesos;
		int *Ix;
		int *Fx;
		
		int HaloPressure;
		int Halo;
		
		int HP;

		// Constructor de la clase
		Parallel(ReadData);
		
		// Functions of the class

			// Initial functions
			void Rango_Procesos();
        	void Total_Procesos();
			void AllocateMemory(Memory);
			void WorkSplit(int, int*, int*);
			void RunParallel(Memory);

			// Communication

				// Local Communication
				void CommunicateDataLP(double*, double*);
				void CommunicateDataLU(double*, double*);
				void CommunicateDataLV(double*, double*);
				void CommunicateDataLW(double*, double*);

				// Global Communication
				void SendMatrixToZeroMP(double*, double*);
				void SendMatrixToZeroMU(double*, double*);
				void SendMatrixToZeroMV(double*, double*);
				void SendMatrixToZeroMW(double*, double*);

};
