//------------------------------------------------------------------------------------------------//
//                     CPP FILE FOR CFD SOLVER CLASS MEMORY ALLOCATION                            //
//------------------------------------------------------------------------------------------------//

// Function to allocate memory for Poisson Coefficients
void CFD_Solver::Allocate_PoissonCoeffsMemory(Memory M1){

    A.aw = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);
    A.ae = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);

    A.as = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);
    A.an = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);

    A.ah = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);
    A.at = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);

    A.ap = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1); 
    A.bp = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);

}

// Function to allocate memory for the components of velocities structures
void CFD_Solver::Allocate_VelocitiesPartMemory(Memory M1, Velocity_Struct &StructName, int NodeX, int NodesY, int NodesZ){

    StructName.Predictor = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);

    StructName.Pres = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);
    StructName.Fut = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);

    StructName.ContributionPast = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);
    StructName.ContributionPres = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);

    StructName.Convective = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);
    StructName.Diffusive = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);

    StructName.Boussinesq = M1.AllocateDouble(NodeX, NodesY, NodesZ, 1);
}

// Function to allocate memory for the boundary conditions of the velocities
void CFD_Solver::Allocate_VelocitiesBoundaryConditionsMemory(Memory M1, Velocity_Struct &StructName, int NodesX, int NodesY, int NodesZ){

    if (Rango == 0){
		StructName.Left = M1.AllocateDouble(1, NodesY, NodesZ, 1);
    }
    else if (Rango == Procesos - 1){
        StructName.Right = M1.AllocateDouble(1, NodesY, NodesZ, 1);
    }

    StructName.Bottom = M1.AllocateDouble(NodesX, 1, NodesZ, 1);
    StructName.Top = M1.AllocateDouble(NodesX, 1, NodesZ, 1);

    StructName.Here = M1.AllocateDouble(NodesX, NodesY, 1, 1);
    StructName.There = M1.AllocateDouble(NodesX, NodesY, 1, 1);

}

// Function to allocate all the memory needed for CFD_Solver class
void CFD_Solver::Allocate_VelocitiesMemory(Memory M1){

    // Velocity U
    Allocate_VelocitiesPartMemory(M1, U, Fx[Rango] - Ix[Rango] + 2*Halo + 1, NY + 2*Halo, NZ + 2*Halo);
    Allocate_VelocitiesBoundaryConditionsMemory(M1, U, Fx[Rango] - Ix[Rango] + 2, NY, NZ);

    // Velocity V
    Allocate_VelocitiesPartMemory(M1, V, Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo + 1, NZ + 2*Halo);
    Allocate_VelocitiesBoundaryConditionsMemory(M1, V, Fx[Rango] - Ix[Rango] + 2, NY, NZ);

    // Velocity W
    Allocate_VelocitiesPartMemory(M1, W, Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, NZ + 2*Halo + 1);
    Allocate_VelocitiesBoundaryConditionsMemory(M1, W, Fx[Rango] - Ix[Rango] + 2, NY, NZ);
  
}

// Function to allocate memory for the pressure matrix
void CFD_Solver::Allocate_PressureMemory(Memory M1){
    P.Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    P.Sup = M1.AllocateDouble(Fx[Rango] - Ix[Rango], NY, NZ, 1);
}

// Function to allocate memory for the energy equation resolution
void CFD_Solver::Allocate_EnergyMemory(Memory M1){

    T.Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    T.Fut = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);

    T.ContributionPast = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    T.ContributionPres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);

    T.Convective = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    T.Diffusive = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    
    if (Rango == 0){
		T.Left = M1.AllocateDouble(1, NY, NZ, 1);
    }
    else if (Rango == Procesos - 1){
        T.Right = M1.AllocateDouble(1, NY, NZ, 1);
    }

    T.Bottom = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, 1, NZ, 1);
    T.Top = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, 1, NZ, 1);

    T.Here = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, NY, 1, 1);
    T.There = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, NY, 1, 1);

}

// Function to allocate memory for the global matrix in core 0
void CFD_Solver::Allocate_GlobalMemory(Memory M1){

    Global.P = M1.AllocateDouble(NX + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    
    Global.U = M1.AllocateDouble(NX + 2*HP + 1, NY + 2*Halo, NZ + 2*Halo, 1);
    Global.V = M1.AllocateDouble(NX + 2*HP, NY + 2*Halo + 1, NZ + 2*Halo, 1);
    Global.W = M1.AllocateDouble(NX + 2*HP, NY + 2*Halo, NZ + 2*Halo + 1, 1);

    Global.T = M1.AllocateDouble(NX + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
}

// Function to delete all the memory of a velocity structure
void CFD_Solver::Delete_VelocityMemory(Velocity_Struct &StructName){

    delete [] StructName.Predictor;

    delete [] StructName.Pres;
    delete [] StructName.Fut;

    delete [] StructName.ContributionPast;
    delete [] StructName.ContributionPres;

    delete [] StructName.Convective;
    delete [] StructName.Diffusive;

    delete [] StructName.Bottom;
    delete [] StructName.Top;

    delete [] StructName.Here;
    delete [] StructName.There;

    delete [] StructName.Left;
    delete [] StructName.Right;

}