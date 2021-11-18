//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR PARALLEL PROGRAMMING CLASS                            //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"

using namespace std;

#define LM(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo) * ((i) - Ix + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo) + ((Fx - Ix + 2*Halo) * (NY + 2*Halo) * (NZ + 2*Halo)) * (dim)

Parallel::Parallel(){
	
    // Data necessary
    
        // Meshing data
        NX = R1.ProblemNumericalData[1];
        NY = R1.ProblemNumericalData[2];
        NZ = R1.ProblemNumericalData[3];

        // Parallel computing data
		Halo = 2;

        

}

void Parallel::Rango_Procesos(){
	int a;
	a = MPI_Comm_rank(MPI_COMM_WORLD, &Rango);
}

void Parallel::Total_Procesos(){
	int a;
	a = MPI_Comm_size(MPI_COMM_WORLD, &Procesos);
}

void Parallel::AllocateMemory(Memory M1){
    Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);
}

void Parallel::WorkSplit(int NX, int &Ix, int &Fx){
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
    MPI_Barrier(MPI_Comm communicator);
    for (i = 0; i < Procesos; i++){
        if(i != Rango){
            MPI_Recv(Ix[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ST);
        }
    }

    MPI_Barrier(MPI_Comm communicator);

    // Send Fx values to every processor
    for (i = 0; i < Procesos; i++){
        if (i != Rango){
            MPI_Send(&Fx[Rango], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_Comm communicator);
    for (i = 0; i < Procesos; i++){
        if(i != Rango){
            MPI_Recv(Fx[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ST);
        }
    }

}

// Local Matrix Communication
void Parallel::CommunicateLocalMatrix(double *LocalSend, double *LocalReceive, int Ix, int Fx){
MPI_Status ST;	

    // Send Everything to the right
	if(Rango != Procesos - 1){
		MPI_Send(&LocalSend[LM(Fx - Halo, 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD);
	}

	if(Rango != 0){
		MPI_Recv(&LocalReceive[LM(Ix - Halo, 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rank-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

    // Send Everything to the left
	if(Rango != 0){
		MPI_Send(&LocalSend[LM(Ix, 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD);
	}

	if(Rango != Procesos - 1){
		MPI_Recv(&LocalReceive[LM(Fx, 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rank+1, 0, MPI_COMM_WORLD, &ST);
	}
	
}