//------------------------------------------------------------------------------------------------//
//                     CPP FILE FOR SPECIES SOLVER CLASS MEMORY ALLOCATION                        //
//------------------------------------------------------------------------------------------------//

// Function to allocate memory for each specie on the solver
void Species_Solver::Allocate_StructSpecies(Memory M1, int SP){

    // Mass Fraction fields maps in each time step
    Species[SP].C_Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 1);
    Species[SP].C_Fut = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 1);

    // Equation terms
    Species[SP].Convective = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);
    Species[SP].Diffusive = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);

    Species[SP].ContributionPres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);
    Species[SP].ContributionPast = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);

    Species[SP].D_ab = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);

    if (Rango == 0){
		Species[SP].Left = M1.AllocateDouble(1, NY, NZ, 1);
    }
    else if (Rango == Procesos - 1){
        Species[SP].Right = M1.AllocateDouble(1, NY, NZ, 1);
    }

    Species[SP].Bottom = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, 1, NZ, 1);
    Species[SP].Top = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, 1, NZ, 1);

    Species[SP].Here = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, NY, 1, 1);
    Species[SP].There = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, NY, 1, 1);
          
    if (Rango == 0){
      Species[SP].Global = M1.AllocateDouble(NX + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    }

}

