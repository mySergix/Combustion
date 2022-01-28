//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR PARALLEL PROGRAMMING CLASS                            //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"



using namespace std;

Parallel::Parallel(ReadData R1){

		NX = R1.ProblemNumericalData[2];
	    NY = R1.ProblemNumericalData[3];
	    NZ = R1.ProblemNumericalData[4];

		Halo = 2;
		HP = 2;
}

// Additional files of the class
#include "Matrix_Index.cpp"
#include "Parallel_Communication.cpp"

// Function to get the rank of each core
void Parallel::Rango_Procesos(){
	int a;
	a = MPI_Comm_rank(MPI_COMM_WORLD, &Rango);
}

// Function to get the total number of cores
void Parallel::Total_Procesos(){
	int a;
	a = MPI_Comm_size(MPI_COMM_WORLD, &Procesos);
}

// Function to allocate memory for all cores Ix and Fx
void Parallel::AllocateMemory(Memory M1){
    Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);
}

// Function to perform the worksplit between all cores
void Parallel::WorkSplit(int NX, int *Ix, int *Fx){
MPI_Status ST;
int i;
int Intervalo, Residuo;
int p;

	Intervalo = NX/Procesos;
	Residuo = NX%Procesos;

	if(Rango != Procesos-1){
		Ix[Rango] = Rango*Intervalo;
		Fx[Rango] = (Rango+1)*Intervalo;
	}
	else{
		Ix[Rango] = Rango*Intervalo;
		Fx[Rango] = (Rango+1)*Intervalo + Residuo;
	}

    // Send Ix values to every processor
    for (i = 0; i < Procesos; i++){
        if (i != Rango){
            MPI_Send(&Ix[Rango], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < Procesos; i++){
        if(i != Rango){
            MPI_Recv(&Ix[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ST);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Send Fx values to every processor
    for (i = 0; i < Procesos; i++){
        if (i != Rango){
            MPI_Send(&Fx[Rango], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < Procesos; i++){
        if(i != Rango){
            MPI_Recv(&Fx[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ST);
        }
    }

}

// Function to execute the parallel class
void Parallel::RunParallel(Memory M1){

    Rango_Procesos();
    Total_Procesos();
    AllocateMemory(M1);
    WorkSplit(NX, Ix, Fx);
    
}