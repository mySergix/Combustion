//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS UTILITIES                            //
//------------------------------------------------------------------------------------------------//

// Function to set the initial fields
void CFD_Solver::Get_InitialConditions(){
int i, j, k;

	// Pressure
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				P.Pres[LP(i,j,k,0)] = 0.0;
			}
		}
	}
	
	// Velocity U
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){	
				U.Pres[LU(i,j,k,0)] = 0.0;
				U.Fut[LU(i,j,k,0)] = 0.0;

				U.Convective[LU(i,j,k,0)] = 0.0;
				U.Diffusive[LU(i,j,k,0)] = 0.0;	
			}
		}
	}
	
    // Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){	
				V.Pres[LV(i,j,k,0)] = 0.0;
				V.Fut[LV(i,j,k,0)] = 0.0;

				V.Convective[LV(i,j,k,0)] = 0.0;
				V.Diffusive[LV(i,j,k,0)] = 0.0;	
			}
		}
	}

    // Velocity W
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ + 1; k++){	
				W.Pres[LW(i,j,k,0)] = 0.0;
				W.Fut[LW(i,j,k,0)] = 0.0;

				W.Convective[LW(i,j,k,0)] = 0.0;
				W.Diffusive[LW(i,j,k,0)] = 0.0;	
			}
		}
	}
	
}

//Cálculo del tiempo entre Steps
void CFD_Solver::Get_StepTime(Mesher MESH, Parallel P1){
int i, j , k;
DeltaT = 1000.0;
double Tpar = 0.20;

	//Comparacion 
	//a*(a>b) + b*(b>a)
	//b+=(a-b)*(a>b)

	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
	
				// CFL Velocities (U, V, W)
	
				// Velocity U
                DeltaT += (abs(MESH.DeltasMP[LP(i,j,k,0)] / (0.50 * (U.Pres[LU(i,j,k,0)] + U.Pres[LU(i+1,j,k,0)]) + 1e-10)) - DeltaT) * (abs(MESH.DeltasMP[LP(i,j,k,0)] / (0.50 * (U.Pres[LU(i,j,k,0)] + U.Pres[LU(i+1,j,k,0)]) + 1e-10)) <= DeltaT);

				// Velocity V
				DeltaT += (abs(MESH.DeltasMP[LP(i,j,k,1)] / (0.50 * (V.Pres[LV(i,j,k,0)] + V.Pres[LV(i,j+1,k,0)]) + 1e-10)) - DeltaT) * (abs(MESH.DeltasMP[LP(i,j,k,1)] / (0.50 * (V.Pres[LV(i,j,k,0)] + V.Pres[LV(i,j+1,k,0)]) + 1e-10)) <= DeltaT);

				//Velocidad W
				DeltaT += (abs(MESH.DeltasMP[LP(i,j,k,2)] / (0.50 * (W.Pres[LW(i,j,k,0)] + W.Pres[LW(i,j,k+1,0)]) + 1e-10)) - DeltaT) * (abs(MESH.DeltasMP[LP(i,j,k,2)] / (0.50 * (W.Pres[LW(i,j,k,0)] + W.Pres[LW(i,j,k+1,0)]) + 1e-10)) <= DeltaT);

			
				//CFL Difusivo Velocidades
                
				//Difusivo U
				DeltaT += ((Tpar*Rho*pow(MESH.DeltasMP[LP(i,j,k,0)],2.0))/(mu + 1e-10) - DeltaT) * ((Tpar*Rho*pow(MESH.DeltasMP[LP(i,j,k,0)],2.0))/(mu + 1e-10) <= DeltaT);

				//Difusivo V
				DeltaT += ((Tpar*Rho*pow(MESH.DeltasMP[LP(i,j,k,1)],2.0))/(mu + 1e-10) - DeltaT) * ((Tpar*Rho*pow(MESH.DeltasMP[LP(i,j,k,1)],2.0))/(mu + 1e-10) <= DeltaT);

				//Difusivo W
				DeltaT += ((Tpar*Rho*pow(MESH.DeltasMP[LP(i,j,k,2)],2.0))/(mu + 1e-10) - DeltaT) * ((Tpar*Rho*pow(MESH.DeltasMP[LP(i,j,k,2)],2.0))/(mu + 1e-10) <= DeltaT);

			}	
		}
	}
	
	MPI_Allreduce(&DeltaT, &DeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

// Function to calculate the convective scheme
inline double CFD_Solver::ConvectiveScheme(double CoordObjetivo, double Velocity, double Coord1, double Phi1, double Coord2, double Phi2, double Coord3, double Phi3, double Coord4, double Phi4, string Esquema){

double PhiObjetivo;

double CoordD;
double PhiD;

double CoordC;
double PhiC;

double CoordU;
double PhiU;

	if (Velocity <= 0.0 || (Phi1 == 0.0 && Coord1 == 0.0)){

		CoordD = Coord2;
		PhiD = Phi2;
		CoordC = Coord3;
		PhiC = Phi3;
		CoordU = Coord4;
		PhiU = Phi4;

	}
	else if(Velocity > 0.0 || (Phi4 == 0.0 && Coord4 == 0.0)){

		CoordD = Coord3;
		PhiD = Phi3;
		CoordC = Coord2;
		PhiC = Phi2;
		CoordU = Coord1;
		PhiU = Phi1;

	}

	//Adimensionalizacion
	double PhiAdimC = (PhiC - PhiU)/(PhiD - PhiU);
	double AdimCoordC = (CoordC - CoordU)/(CoordD - CoordU);
	double AdimCoordE = (CoordObjetivo - CoordU)/(CoordD - CoordU);

	//Evaluacion
	double PhiF;

	if (PhiD == PhiU){
		PhiObjetivo = PhiD;
	}
	else{

        // CDS
        PhiF = ((AdimCoordE - AdimCoordC) / (1 - AdimCoordC)) + ((AdimCoordE - 1) / (AdimCoordC - 1)) * PhiAdimC;

        // QUICK
		//PhiF = AdimCoordE + (((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0))))*(PhiAdimC - AdimCoordC);

		//Dimensionalizacion
		PhiObjetivo = PhiU + (PhiD - PhiU)*PhiF;
	}

	return PhiObjetivo;

}

// Function to calculate the velocity contributions and predictor velocities
void CFD_Solver::Get_ContributionsPredictors(){
int i, j, k;

    // Velocity U
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				U.ContributionPres[LU(i,j,k,0)] = U.Diffusive[LU(i,j,k,0)] - U.Convective[LU(i,j,k,0)];
			}
		}
	}

	// Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY + 1; j++){
		    for(k = 0; k < NZ; k++){
				V.ContributionPres[LV(i,j,k,0)] = V.Diffusive[LV(i,j,k,0)] - V.Convective[LV(i,j,k,0)];
			}
		}
	}

    // Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ + 1; k++){
				W.ContributionPres[LW(i,j,k,0)] = W.Diffusive[LW(i,j,k,0)] - W.Convective[LW(i,j,k,0)];
			}
		}
	}

	// Velocity Predictor U
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				U.Predictor[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + DeltaT*(1.50*U.ContributionPres[LU(i,j,k,0)] - 0.50*U.ContributionPast[LU(i,j,k,0)]);
			}
		}
	}

    // Velocity Predictor V
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY + 1; j++){
		    for(k = 0; k < NZ; k++){
				V.Predictor[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + DeltaT*(1.50*V.ContributionPres[LV(i,j,k,0)] - 0.50*V.ContributionPast[LV(i,j,k,0)]);
			}
		}
	}

    // Velocity Predictor W
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ + 1; k++){
				W.Predictor[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + DeltaT*(1.50*W.ContributionPres[LW(i,j,k,0)] - 0.50*W.ContributionPast[LW(i,j,k,0)]);
			}
		}
	}

}

// Function to calculate the divergence of the predictors velocities (bp)
void CFD_Solver::Get_PredictorsDivergence(Mesher MESH){
int i, j, k;

    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				A.bp[LA(i,j,k,0)] = - (Rho/DeltaT)*(
								    + MESH.SupMP[LP(i,j,k,0)] * (U.Predictor[LU(i+1,j,k,0)] - U.Predictor[LU(i,j,k,0)])
								    + MESH.SupMP[LP(i,j,k,1)] * (V.Predictor[LV(i,j+1,k,0)] - V.Predictor[LV(i,j,k,0)])
								    + MESH.SupMP[LP(i,j,k,2)] * (W.Predictor[LW(i,j,k+1,0)] - W.Predictor[LW(i,j,k,0)])
								    );
			}
		}
	}

}

// Function to Check the difference between time steps
void CFD_Solver::Get_Stop(){
int i, j, k;
MaxDiffGlobal = 0.0;

	// Velocity U
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				MaxDiffGlobal += (abs((U.Fut[LU(i,j,k,0)] - U.Pres[LU(i,j,k,0)])/(U.Pres[LU(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((U.Fut[LU(i,j,k,0)] - U.Pres[LU(i,j,k,0)])/(U.Pres[LU(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	// Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY + 1; j++){
		    for(k = 0; k < NZ; k++){
				MaxDiffGlobal += (abs((V.Fut[LV(i,j,k,0)] - V.Pres[LV(i,j,k,0)])/(V.Pres[LV(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((V.Fut[LV(i,j,k,0)] - V.Pres[LV(i,j,k,0)])/(V.Pres[LV(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	// Velocity W
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ + 1; k++){
				MaxDiffGlobal += (abs((W.Fut[LW(i,j,k,0)] - W.Pres[LW(i,j,k,0)])/(W.Pres[LW(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((W.Fut[LW(i,j,k,0)] - W.Pres[LW(i,j,k,0)])/(W.Pres[LW(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	MPI_Allreduce(&MaxDiffGlobal, &MaxDiffGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

}

// Function to update the matrix fields
void CFD_Solver::Get_Update(){
int i, j, k;

	// Velocity U
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				U.Pres[LU(i,j,k,0)] = U.Fut[LU(i,j,k,0)];
                U.ContributionPast[LU(i,j,k,0)] = U.ContributionPres[LU(i,j,k,0)];
            }
		}
	}

	// Velocity V
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY + 1; j++){
		    for(k = 0; k < NZ; k++){
				V.Pres[LV(i,j,k,0)] = V.Fut[LV(i,j,k,0)];
                V.ContributionPast[LV(i,j,k,0)] = V.ContributionPres[LV(i,j,k,0)];
            }
		}
	}

	// Velocity W
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ + 1; k++){
				W.Pres[LW(i,j,k,0)] = W.Fut[LW(i,j,k,0)];
                W.ContributionPast[LW(i,j,k,0)] = W.ContributionPres[LW(i,j,k,0)];
            }
		}
	}

}

// Function to calculate the new corrected velocities
void CFD_Solver::Get_Velocities(Mesher MESH, Parallel P1){
int i, j, k;

	// Calculated Pressure Halo Communication
	P1.CommunicateDataLP(P.Pres, P.Pres);

	// Velocity U
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U.Fut[LU(i,j,k,0)] = U.Predictor[LU(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]);
			}
		}
	}

    if (Rango == 0){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U.Fut[LU(0,j,k,0)] = U.Left[LEFT(0,j,k)];
			}
		}
    }
    else if (Rango == Procesos - 1){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U.Fut[LU(NX,j,k,0)] = U.Right[RIGHT(NX,j,k)];
			}
		}
    }

    // Velocity V
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(k = 0; k < NZ; k++){

            // Bottom
            V.Fut[LV(i,0,k,0)] = V.Bottom[BOTTOM(i,0,k)];

            // Top
            V.Fut[LV(i,NY,k,0)] = V.Top[TOP(i,NY,k)];

            for(j = 1; j < NY; j++){ 		
				V.Fut[LV(i,j,k,0)] = V.Predictor[LV(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]);
			}
		}
	}

    // Velocity W
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){

            // Here
            W.Fut[LW(i,j,0,0)] = W.Here[HERE(i,j,0)];

            // There
            W.Fut[LW(i,j,NZ,0)] = W.There[THERE(i,j,NZ)];

		    for(k = 1; k < NZ; k++){		
				W.Fut[LW(i,j,k,0)] = W.Predictor[LW(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,0)]);
			}
		}
	}

}
