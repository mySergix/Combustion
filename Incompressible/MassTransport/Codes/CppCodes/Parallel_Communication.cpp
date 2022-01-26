//------------------------------------------------------------------------------------------------//
//                      CPP FILE FOR PARALLEL COMMUNICATION FUNCTIONS                             //
//------------------------------------------------------------------------------------------------//

// Function to communicate collocated matrix
void Parallel::CommunicateDataLP(double *LocalSend, double *LocalReceive){
MPI_Status ST;	

	if(Rango != Procesos - 1){
		MPI_Send(&LocalSend[LP(Fx[Rango] - HP, 0, 0, 0)], (HP)*(NY)*(NZ), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD);
	}

	if(Rango != 0){
		MPI_Recv(&LocalReceive[LP(Ix[Rango] - HP, 0, 0, 0)], (HP)*(NY)*(NZ), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rango != 0){
		MPI_Send(&LocalSend[LP(Ix[Rango], 0, 0, 0)], (HP)*(NY)*(NZ), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD);
	}

	if(Rango != Procesos - 1){
		MPI_Recv(&LocalReceive[LP(Fx[Rango], 0, 0, 0)], (HP)*(NY)*(NZ), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD, &ST);
	}
	
}

// Function to communicate staggered matrix U
void Parallel::CommunicateDataLU(double *LocalSend, double *LocalReceive){
MPI_Status ST;

	if(Rango != Procesos - 1){
		MPI_Send(&LocalSend[LU(Fx[Rango] - Halo, - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD);
	}

	if(Rango != 0){
		MPI_Recv(&LocalReceive[LU(Ix[Rango] - Halo, - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rango != 0){
		MPI_Send(&LocalSend[LU(Ix[Rango] + 1, - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD);
	}

	if(Rango != Procesos - 1){
		MPI_Recv(&LocalReceive[LU(Fx[Rango] + 1, - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD, &ST);
	}
	
}

// Function to communicate staggered matrix V
void Parallel::CommunicateDataLV(double *LocalSend, double *LocalReceive){
MPI_Status ST;

	if(Rango != Procesos - 1){
		MPI_Send(&LocalSend[LV(Fx[Rango] - Halo, - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo + 1)*(NZ + 2*Halo), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD);
	}

	if(Rango != 0){
		MPI_Recv(&LocalReceive[LV(Ix[Rango] - Halo, - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo + 1)*(NZ + 2*Halo), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rango != 0){
		MPI_Send(&LocalSend[LV(Ix[Rango], - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo + 1)*(NZ + 2*Halo), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD);
	}

	if(Rango != Procesos - 1){
		MPI_Recv(&LocalReceive[LV(Fx[Rango], - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo + 1)*(NZ + 2*Halo), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD, &ST);
	}

}

// Function to communicate staggered matrix W
void Parallel::CommunicateDataLW(double *LocalSend, double *LocalReceive){
MPI_Status ST;
	
	if(Rango != Procesos - 1){
		MPI_Send(&LocalSend[LW(Fx[Rango] - Halo,- Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo + 1), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD);
	}

	if(Rango != 0){
		MPI_Recv(&LocalReceive[LW(Ix[Rango] - Halo, - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo + 1), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD, &ST);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(Rango != 0){
		MPI_Send(&LocalSend[LW(Ix[Rango], - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo + 1), MPI_DOUBLE, Rango-1, 0, MPI_COMM_WORLD);
	}

	if(Rango != Procesos - 1){
		MPI_Recv(&LocalReceive[LW(Fx[Rango], - Halo, - Halo, 0)], (Halo)*(NY + 2*Halo)*(NZ + 2*Halo + 1), MPI_DOUBLE, Rango+1, 0, MPI_COMM_WORLD, &ST);
	}

}

// Function to communicate Collocated P Matrix to Processor 0
void Parallel::SendMatrixToZeroMP(double *LocalMatrix, double *GlobalMatrix){
int i, j, k, p;
MPI_Status ST;

	if(Rango != 0){
		MPI_Send(&LocalMatrix[LP(Ix[Rango], - HP, - HP, 0)], (Fx[Rango] - Ix[Rango])*(NY + 2*HP)*(NZ + 2*HP), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rango == 0){
		for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					GlobalMatrix[GP(i,j,k,0)] = LocalMatrix[LP(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			MPI_Recv(&GlobalMatrix[GP(Ix[p], - HP, - HP, 0)], (Fx[p] - Ix[p])*(NY + 2*HP)*(NZ+ 2*HP), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);	

}

// Function to communicate Staggered U Matrix to Processor 0
void Parallel::SendMatrixToZeroMU(double *LocalMatrix, double *GlobalMatrix){
int i, j, k, p;
MPI_Status ST;

	if(Rango != 0){
		MPI_Send(&LocalMatrix[LU(Ix[Rango], - Halo, - Halo, 0)], (Fx[Rango] - Ix[Rango] + 1)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rango == 0){
		for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
			for(j = - Halo; j < NY + Halo; j++){
				for(k = - Halo; k < NZ + Halo; k++){
					GlobalMatrix[GU(i,j,k,0)] = LocalMatrix[LU(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			MPI_Recv(&GlobalMatrix[GU(Ix[p], - Halo, - Halo, 0)], (Fx[p] - Ix[p] + 1)*(NY + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}
	
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

// Function to communicate Staggered U Matrix to Processor 0
void Parallel::SendMatrixToZeroMV(double *LocalMatrix, double *GlobalMatrix){
int i, j, k, p;
MPI_Status ST;

	if(Rango != 0){
		MPI_Send(&LocalMatrix[LV(Ix[Rango], - Halo, - Halo, 0)], (Fx[Rango] - Ix[Rango])*(NY + 1 + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rango == 0){
		for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = - Halo; j < NY + 1 + Halo; j++){
				for(k = - Halo; k < NZ + Halo; k++){
					GlobalMatrix[GV(i,j,k,0)] = LocalMatrix[LV(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){	
			MPI_Recv(&GlobalMatrix[GV(Ix[p], - Halo, - Halo, 0)], (Fx[p] - Ix[p])*(NY + 1 + 2*Halo)*(NZ + 2*Halo), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);
		}
	
	}

	MPI_Barrier(MPI_COMM_WORLD);

}

// Function to communicate Staggered U Matrix to Processor 0
void Parallel::SendMatrixToZeroMW(double *LocalMatrix, double *GlobalMatrix){
int i, j, k, p;
MPI_Status ST;

	if(Rango != 0){
		MPI_Send(&LocalMatrix[LW(Ix[Rango], - Halo, - Halo, 0)], (Fx[Rango] - Ix[Rango])*(NY + 2*Halo)*(NZ + 1 + 2*Halo), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(Rango == 0){
		
		for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = - Halo; j < NY + Halo; j++){
				for(k = - Halo; k < NZ + 1 + Halo; k++){
					GlobalMatrix[GW(i,j,k,0)] = LocalMatrix[LW(i,j,k,0)];
				}		
			}
		}
		
		for(p = 1; p < Procesos; p++){
			MPI_Recv(&GlobalMatrix[GW(Ix[p], - Halo, - Halo, 0)], (Fx[p] - Ix[p])*(NY + 2*Halo)*(NZ + 1 + 2*Halo), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &ST);	
		}
	
	}

	MPI_Barrier(MPI_COMM_WORLD);

}