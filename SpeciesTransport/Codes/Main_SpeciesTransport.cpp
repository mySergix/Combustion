#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include </usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h>



#include "HeaderCodes/Memory.h"
#include "HeaderCodes/ReadData.h"
//#include "HeaderCodes/CanteraCombustion.h"

// Clase para leer los input de Cantera
// Clase para la programacion en paralelo
// Clase para el mallador
// Clase para el solver
// Clase para el postproceso




int main(int argc, char* argv[]){

MPI_Init(&argc, &argv);
cout << endl;

Memory M1;

ReadData R1(M1);
R1.ReadInputs(M1);


MPI_Finalize();

return 0;

}