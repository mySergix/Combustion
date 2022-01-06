//------------------------------------------------------------------------------------------------//
//                     CPP FILE FOR CFD SOLVER CLASS MEMORY ALLOCATION                            //
//------------------------------------------------------------------------------------------------//

// Function to allocate memory for the maps fields and wall values
void CFD_Solver::Allocate_Struct_MapFields(Memory M1, Property_Struct &PropertyName){

    // Fields maps in each time step
    PropertyName.Past = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Fut = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    // Properties values at the walls
    PropertyName.Wall_U = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1); 
    PropertyName.Wall_V = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Wall_W = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

}

// Function to allocate memory for the velocity gradients
void CFD_Solver::Allocate_Struct_VelocityGradients(Memory M1, Property_Struct &PropertyName){

    // Properties gradients at the walls
    PropertyName.Wall_U_Diff_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Wall_U_Diff_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Wall_U_Diff_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    PropertyName.Wall_V_Diff_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Wall_V_Diff_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Wall_V_Diff_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);

    PropertyName.Wall_W_Diff_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);
    PropertyName.Wall_W_Diff_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);
    PropertyName.Wall_W_Diff_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

}

// Function to allocate memory for the boundary conditions of a property
void CFD_Solver::Allocate_Struct_BoundaryConditions(Memory M1, Property_Struct &PropertyName){

    // Properties Boundary Conditions 
    PropertyName.Bottom = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, 1, NZ + 2*Halo, 1);
    PropertyName.Top = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, 1, NZ + 2*Halo, 1);

    PropertyName.Here = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, 1, 1);
    PropertyName.There = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, 1, 1);

    if (Rango == 0){
        PropertyName.Left = M1.AllocateDouble(1, NY + 2*Halo, NZ + 2*Halo, 1);     
    }
    if (Rango == Procesos - 1){
        PropertyName.Right = M1.AllocateDouble(1, NY + 2*Halo, NZ + 2*Halo, 1);  
    }

}

// Function to allocate the equations contributions
void CFD_Solver::Allocate_Struct_Contributions(Memory M1, Property_Struct &PropertyName){

    // Properties equations contributions
    PropertyName.ContributionPast = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.ContributionPres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

// Function to allocate memory for velocity equations terms
void CFD_Solver::Allocate_Struct_VelocityEqsTerms(Memory M1, Property_Struct &PropertyName){

    // Velocities equations terms
    PropertyName.Convective = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    PropertyName.Diffusive = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

// Function to allocate memory for pressure structure
void CFD_Solver::Allocate_Struct_Pressure(Memory M1){

    // Pressure matrix
    Pressure.Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Pressure.Gradient_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Pressure.Gradient_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Pressure.Gradient_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

// Function to allocate memory for the viscous stresses
void CFD_Solver::Allocate_Struct_Stresses(Memory M1){

    // Viscous stresses terms
    Stress.mu_Visc = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    Stress.Divergence = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    Stress.Tau_xx_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_yy_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_zz_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    Stress.Tau_yx_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_zx_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    Stress.Tau_xy_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_zy_Z = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    Stress.Tau_yz_Y = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Stress.Tau_xz_X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

// Function to allocate memory for aditional energy equation terms
void CFD_Solver::Allocate_Struct_EnergyEqTerms(Memory M1){

    // Energy equation aditional terms
    Lambda = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    FourierDiffusion = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1); 
    EnthalpyDiffusion = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

// Function to allocate memory of global map fields
void CFD_Solver::Allocate_Struct_Global(Memory M1){

    // Global map fields
    GlobalMatrix.Density = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    GlobalMatrix.U = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    GlobalMatrix.V = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    GlobalMatrix.W = M1.AllocateDouble(NX + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

}

