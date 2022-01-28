//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the solver of the simulation
void CFD_Solver::RunSolver(Memory M1, Parallel P1, Mesher MESH, PostProcessing POST1){	
int i, j, k;

	double Time = 0.0;
	int Step = 0;

	char FileName_1[300]; 

	MaxDiffGlobal = 2.0*ConvergenciaGlobal;

    // Memory Allocation
	Allocate_VelocitiesMemory(M1);
    Allocate_PressureMemory(M1);
    Allocate_PoissonCoeffsMemory(M1);
	Allocate_EnergyMemory(M1);
	if (Rango == 0){ Allocate_GlobalMemory(M1); }

    // Initial settings and calculations
    Get_InitialConditions();
    Get_PoissonCoefficients(MESH);
	PrintTxt();
	
    Get_InitialBoundaryConditions();
	Get_StaticHalos();

	while(MaxDiffGlobal >= ConvergenciaGlobal){
		
		Step++;
		
		Get_UpdateBoundaryConditions();
		Get_StaticHalos();
        Get_UpdateHalos();

		P1.CommunicateDataLU(U.Pres, U.Pres);
		P1.CommunicateDataLV(V.Pres, V.Pres);
		P1.CommunicateDataLW(W.Pres, W.Pres);
		
		P1.CommunicateDataLP(T.Pres, T.Pres);

		Get_StepTime(MESH, P1); // Time Step Calculation
		Time += DeltaT;
        
        Get_DiffusionU(MESH);
        Get_DiffusionV(MESH);
        Get_DiffusionW(MESH);
        
		Get_ConvectionU(MESH);
		Get_ConvectionV(MESH);
		Get_ConvectionW(MESH);

		if (Problema == 2){
			Get_BoussinesqV(MESH);
		}

        Get_ContributionsPredictors();
        Get_PredictorsDivergence(MESH);

        Get_GaussSeidel(P1);
		Get_Velocities(MESH, P1);

		Get_DiffusionEnergy(MESH);
        //Get_ConvectionEnergy(MESH);
        
        Get_Temperature();

		if(Step%100 == 0){

			// Communication to global matrix
			P1.SendMatrixToZeroMP(P.Pres, Global.P);
			P1.SendMatrixToZeroMP(T.Fut, Global.T);
			P1.SendMatrixToZeroMU(U.Fut, Global.U);
			P1.SendMatrixToZeroMV(V.Fut, Global.V);
			P1.SendMatrixToZeroMW(W.Fut, Global.W);

			Get_Stop();
			if(Rango == 0){
				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);
				if (Problema == 2){
					sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
					POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
				}
				sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
				POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<endl;
			}

		}

		Get_Update();
		
	}
	
	if(Rango == 0){
		cout<<"Solver Completed"<<endl;
	}

}
