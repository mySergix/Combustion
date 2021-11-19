//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR PARALLEL PROGRAMMING CLASS                            //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"

using namespace std;

#define LM(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo) + ((Fx[Rango] - Ix[Rango] + 2*Halo) * (NY + 2*Halo) * (NZ + 2*Halo)) * (dim)
#define GM(i,j,k,dim) (NY + 2*Halo) * (NZ + 2*Halo) * ((i) + Halo) + (NZ + 2*Halo) * ((j) + Halo) + ((k) + Halo) + (NY + 2*Halo) * (NZ + 2*Halo) * (NX + 2*Halo) * (dim)

Parallel::Parallel(ReadData R1){
	
    // Data necessary

        // Meshing data
        NX = R1.ProblemNumericalData[0];
        NY = R1.ProblemNumericalData[1];
        NZ = R1.ProblemNumericalData[2];

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

void Parallel::WorkSplit(int NX, int* Ix, int* Fx){
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

// Local Matrix Communication
void Parallel::CommunicateLocalMatrix(double* LocalSend, double* LocalReceive){
MPI_Status ST;	

    // Send Everything to the right
	if(Rango != Procesos - 1){
		MPI_Send(&LocalSend[LM(Fx[Halo] - Halo, 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD);
	}
    
	if(Rango != 0){
		MPI_Recv(&LocalReceive[LM(Ix[Rango] - Halo, 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD, &ST);
	}
    
	MPI_Barrier(MPI_COMM_WORLD);

    // Send Everything to the left
	if(Rango != 0){
		MPI_Send(&LocalSend[LM(Ix[Rango], 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD);
	}
    
	if(Rango != Procesos - 1){
		MPI_Recv(&LocalReceive[LM(Fx[Rango], 0, 0, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD, &ST);
	}
	
}

void Parallel::SendMatrixToZero(double *LocalMatrix, double *GlobalMatrix){
int i, j, k;
MPI_Status ST;

	if(Rango != 0){
		MPI_Send(&LocalMatrix[LM(Ix[Rango], - Halo, - Halo, 0)], (Fx[Rango] - Ix[Rango]) * (NY + 2*Halo) * (NZ + 2*Halo), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rango == 0){

		for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
			for(j = - Halo; j < NY + Halo; j++){
				for(k = - Halo; k < NZ + Halo; k++){
					GlobalMatrix[GM(i,j,k,0)] = LocalMatrix[LM(i,j,k,0)];
				}		
			}
		}
		
		for(i = 1; i < Procesos; i++){
			MPI_Recv(&GlobalMatrix[GM(Ix[i], - Halo, - Halo, 0)], (Fx[i] - Ix[i]) * (NY + 2*Halo) * (NZ + 2*Halo), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &ST);
		}
	
	}

	MPI_Barrier(MPI_COMM_WORLD);	

}

void Parallel::RunParallel(Memory M1){

    Rango_Procesos();
    Total_Procesos();
    AllocateMemory(M1);
    WorkSplit(NX, Ix, Fx);

}