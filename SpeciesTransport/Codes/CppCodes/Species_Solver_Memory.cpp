//------------------------------------------------------------------------------------------------//
//                     CPP FILE FOR SPECIES SOLVER CLASS MEMORY ALLOCATION                        //
//------------------------------------------------------------------------------------------------//

// Function to allocate memory for the maps fields and wall values
void Species_Solver::Allocate_Struct_Species(Memory M1, int SP){

    // Mass Fraction fields maps in each time step
    Species[SP].Y_Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Species[SP].Y_Fut = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    // Mass Fraction values at the walls
    Species[SP].Y_Wall_U = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1); 
    Species[SP].Y_Wall_V = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Species[SP].Y_Wall_W = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    // Molar fraction
    Species[SP].X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    // Diffusion coefficient
    Species[SP].D_am = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    // Diffusion velocities
    Species[SP].U_Diff = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 1 + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Species[SP].V_Diff = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 1 + 2*Halo, NZ + 2*Halo, 1);
    Species[SP].W_Diff = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 1 + 2*Halo, 1);

    // Equation terms
    Species[SP].ContributionPres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);
    Species[SP].ContributionPast = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * Halo, NY + 2*Halo, NZ + 2*Halo, 1);

    // JANAF Terms
    Species[SP].Cp_coeff = M1.AllocateDouble(10, 1, 1, 1);
    Species[SP].h_coeff = M1.AllocateDouble(12, 1, 1, 1);
    Species[SP].mu_coeff = M1.AllocateDouble(8, 1, 1, 1);
    Species[SP].lambda_coeff = M1.AllocateDouble(8, 1, 1, 1);

}