//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR MEMORY ALLOCATION CLASS                            //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

class Memory{	

	public:
		//Constructor de la clase
		Memory();
		
		//Metodos de la clase

			//Métodos de alojamiento de memoria
			int *AllocateInt(int, int, int, int);
			double *AllocateDouble(int, int, int, int);
			
};
