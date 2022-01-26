//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the solver of the simulation
void CFD_Solver::RunSolver(Memory M1, Parallel P1, Mesher MESH){	
int i, j, k;

	double Time = 0.0;
	int Step = 0;

	MaxDiffGlobal = 2.0*ConvergenciaGlobal;

    // Memory Allocation
	Allocate_VelocitiesMemory(M1);
    Allocate_PressureMemory(M1);
    Allocate_PoissonCoeffsMemory(M1);

    // Initial settings and calculations
    Get_InitialConditions();
    Get_PoissonCoefficients(MESH);
    Get_InitialBoundaryConditions();
    
	while(MaxDiffGlobal >= ConvergenciaGlobal){
		
		Step++;
		Get_PeriodicBoundaryConditions(); // Updating periodic Boundary Conditions
		Get_StepTime(MESH, P1); // Time Step Calculation
		Time += DeltaT;
        
        Get_DiffusionU(MESH);
        Get_DiffusionV(MESH);
        Get_DiffusionW(MESH);
        
        Get_ContributionsPredictors();
        Get_PredictorsDivergence(MESH);

        Get_GaussSeidel(P1);
		Get_Velocities(MESH, P1);

		if(Step%10 == 0){

			Get_Stop();
			if(Rango == 0){
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<endl;
			}

		}

		Get_Update();
		
	}
	
	if(Rango == 0){
		cout<<"Solver Completed"<<endl;
	}

}
