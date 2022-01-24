//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the CFD Solver
void CFD_Solver::RunSolver(Memory M1, Parallel P1, Mesher MESH, PostProcess PP1, Species_Solver SPE_S1){
int n;

double Time = 0.0; // Time of the simulation
int Step = 0; // Steps of the simulation
MaxDiff = 2.0 * GlobalConvergence;

    // CFD Solver Preprocess

        // Memory allocation

            // Density
            Allocate_Struct_MapFields(M1, Density);
            Allocate_Struct_BoundaryConditions(M1, Density);
            Allocate_Struct_Contributions(M1, Density);

            // Velocities
            Allocate_VelocityMemory(M1, U); // U Velocity
            Allocate_VelocityMemory(M1, V); // U Velocity
            Allocate_VelocityMemory(M1, W); // U Velocity
        
            // Viscous stresses
            Allocate_Struct_Stresses(M1);

            // Pressure
            Allocate_Struct_Pressure(M1);

            // Energy
            Allocate_Struct_MapFields(M1, Hs);
            Allocate_Struct_BoundaryConditions(M1, Hs);
            Allocate_Struct_Contributions(M1, Hs);
            Allocate_Struct_EnergyEqTerms(M1);

        // Initial field settings
        Set_InitialValues();

        // Initial boundary conditions
        // Analizar y poner bien las boundary conditions


    // Species Solver Preprocess

        // Memory allocation
        for (n = 0; n < N_Species; n++){
            SPE_S1.Allocate_Struct_Species(M1, n);
        }

        // Species Data reading
        SPE_S1.Read_SpeciesName("Species_Data.txt");
        SPE_S1.Read_AllSpeciesData();

        // Initial species field settings
        // Calculos y campos iniciales

        // Initial species boundary conditions
        // Analizar y poner bien las boundary conditions


        while (Step < 10){

            Step++;
            Update_BoundaryConditions(MESH);
            Get_TimeStep(MESH);
            Time += DeltaT;

            // CFD Solver
            
                // Walls' Velocities
                Get_WallsValue_Property(P1, MESH, U); // Velocity U
                Get_WallsValue_Property(P1, MESH, U); // Velocity U
                Get_WallsValue_Property(P1, MESH, U); // Velocity U
                
                // Density
                Get_WallsValue_Property(P1, MESH, Density);
                Get_ConvectiveTerm_Density(MESH);
                Get_TemporalIntegration_Property(Density);

        

            // Species Solver


            if (Step % 1 == 0){
            

                if (Rango == 0){
                    cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiff<<endl;
                }

                MPI_Barrier(MPI_COMM_WORLD);

            }

            UpdateField(Density);

        }

}
