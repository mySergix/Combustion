//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the CFD Solver
void CFD_Solver::RunSolver(Memory M1, Parallel P1, Mesher MESH, PostProcess PP1, Species_Solver SPE_S1){

    double Time = 0.0; // Time of the simulation
    int Step = 0; // Steps of the simulation

    MaxDiff = 2.0 * GlobalConvergence;

    char FileName_1[300];


    // CFD Solver Preprocess

        // Memory allocation
        // Initial field settings
        // Initial boundary conditions


    // Species Solver Preprocess

        // Memory allocation
        // Read species data
        // Initial species field settings
        // Initial species boundary conditions

        

    AllocateMemory(M1); // Allocate memory for the matrix
    Set_InitialValues(); // Get the inital values of the properties
    Get_BoundaryConditions(); // Get the fixed boundary conditions
    
    
    while (MaxDiff >= GlobalConvergence){

        Step++;
        Update_BoundaryConditions(MESH);
        Get_TimeStep(MESH);
        Time += DeltaT;

        CommunicateVelocities(P1, U.Pres, V.Pres, W.Pres);
        Get_WallsVelocities(MESH);
        Get_WallsValue_Scalar(P1, MESH.Node_Mesh, Density);

        // Density
        Get_ConvectiveTerm(MESH, Density);
        Get_StepContribution_Density(Density);
        Get_TemporalIntegration(Density);

        // Velocities
        //Get_Pressure();
        //Get_PressureGradient(MESH);

        // Velocity U
        Get_ConvectiveTerm(MESH, U);
        Get_StepContribution_Velocity(U, Pressure.Gradient_X);

        // Velocity V
        Get_ConvectiveTerm(MESH, V);
        Get_StepContribution_Velocity(V, Pressure.Gradient_Y);

        // Velocity W
        Get_ConvectiveTerm(MESH, W);
        Get_StepContribution_Velocity(W, Pressure.Gradient_Z);

        Get_TemporalIntegration(U);
        Get_TemporalIntegration(V);
        Get_TemporalIntegration(W);

        if (Step % 100 == 0){
            Get_ConvergenceCriteria();
            P1.SendMatrixToZero(Density.Pres, GlobalMatrix.Density);
            
            P1.SendMatrixToZero(U.Pres, GlobalMatrix.U);
            P1.SendMatrixToZero(V.Pres, GlobalMatrix.V);
            P1.SendMatrixToZero(W.Pres, GlobalMatrix.W);

            if (Rango == 0){
                cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiff<<endl;
                sprintf(FileName_1, "Density_Field_Step_%d", Step);
                PP1.GlobalEscalarVTK(MESH, "CombustionResults/", "Densidad", FileName_1, GlobalMatrix.Density, 0);

                sprintf(FileName_1, "Velocities_Field_Step_%d", Step);
                PP1.GlobalVectorialVTK(MESH, "CombustionResults/", "Velocities", FileName_1, GlobalMatrix.U, GlobalMatrix.V, GlobalMatrix.W, 0);
            }

            MPI_Barrier(MPI_COMM_WORLD);	
        }
        
        UpdateField(Density);
        UpdateField(U);
        UpdateField(V);
        UpdateField(W);

    }
    
}
