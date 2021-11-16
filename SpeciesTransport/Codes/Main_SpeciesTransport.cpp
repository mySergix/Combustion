#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

using namespace std;

#include "HeaderCodes/Memory.h"
#include "HeaderCodes/ReadData.h"

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

for (int i = 0; i < 3; i++){
    cout<<"MeshingOptionsData: "<<R1.MeshingOptionsData[i]<<", Number: "<<i + 1<<endl;
}

MPI_Finalize();

return 0;

}