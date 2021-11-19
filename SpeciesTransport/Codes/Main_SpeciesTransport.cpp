//------------------------------------------------------------------------------------------------//
//                             MAIN FILE FOR SPECIES TRANSPORT SIMULATION                         //
//------------------------------------------------------------------------------------------------//

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
#include "HeaderCodes/Parallel.h"
#include "HeaderCodes/Mesher.h"

// Clase para el solver
// Clase para el postproceso

int main(int argc, char* argv[]){

MPI_Init(&argc, &argv);
cout << endl;

Memory M1;

ReadData R1(M1);
R1.ReadInputs(M1);

Parallel P1(R1);
P1.RunParallel(M1);

Mesher MESH(R1, P1);
MESH.RunMesher(M1, P1);

MPI_Finalize();

return 0;

}