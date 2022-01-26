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