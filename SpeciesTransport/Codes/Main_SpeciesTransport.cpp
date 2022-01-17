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
#include <mpi.h>

#include "HeaderCodes/Memory.h"
#include "HeaderCodes/ReadData.h"
#include "HeaderCodes/Parallel.h"
#include "HeaderCodes/Mesher.h"
#include "HeaderCodes/PostProcess.h"
#include "HeaderCodes/Species_Solver.h"
#include "HeaderCodes/CFD_Solver.h"

using namespace std;

int main(int argc, char* argv[]){

MPI_Init(&argc, &argv);
cout << endl;

Memory M1;

ReadData R1(M1);
R1.ReadInputs(M1);

Parallel P1(R1);
P1.RunParallel(M1);

Mesher MESH(M1, R1, P1);
MESH.RunMesher(M1);

PostProcess PP1(M1, R1, P1);

Species_Solver SPE_S1(M1, R1, P1);

CFD_Solver CFD_S1(M1, R1, P1, SPE_S1);



//CFD_S1.RunSolver(M1, P1, MESH, PP1);

MPI_Finalize();

return 0;

}