//------------------------------------------------------------------------------------------------//
//                             MAIN FILE FOR MASS TRANSPORT SIMULATION                            //
//------------------------------------------------------------------------------------------------//

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include "HeaderCodes/Memory.h"
#include "HeaderCodes/ReadData.h"
#include "HeaderCodes/Parallel.h"
#include "HeaderCodes/Mesher.h"
//#include "HeaderCodes/PostProcessing.h"
#include "HeaderCodes/CFD_Solver.h"

using namespace std;

int main(int argc, char* argv[]){

int i;

MPI_Init(&argc, &argv);
cout << endl;

Memory M1;

ReadData R1(M1);
R1.ReadInputs();

Parallel P1(R1);
P1.RunParallel(M1);

Mesher MESH(M1, R1, P1);
MESH.ExecuteMesher(M1);

MPI_Barrier(MPI_COMM_WORLD);

//PostProcessing POST1(M1, R1, MESH, DIRECTORIO);

CFD_Solver CFD_S1(M1, R1, P1);

//CFD_S1.RunSolver(M1, P1, MESH);

MPI_Finalize();

return 0;

}
