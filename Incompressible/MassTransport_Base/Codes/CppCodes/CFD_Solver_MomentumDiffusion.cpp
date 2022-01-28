//------------------------------------------------------------------------------------------------//
//                CPP FILE FOR N-S MOMENTUM EQUATION DIFFUSION CALCULATIONS                       //
//------------------------------------------------------------------------------------------------//

// Function to calculate the diffusion term of U Velocity
void CFD_Solver::Get_DiffusionU(Mesher MESH){
int i, j, k;


	
    // Centro
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 1; j < NY - 1; j++){
		    for(k = 1; k < NZ - 1; k++){	
				U.Diffusive[LU(i,j,k,0)] =  (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j+1,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k+1,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
		}
	}
	
	
    // Top
    j = NY - 1;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(k = 1; k < NZ - 1; k++){
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Top[TOP(i,j,k)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k+1,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }
    }
	
	// Bottom
    j = 0;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(k = 1; k < NZ - 1; k++){
            U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j+1,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k+1,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }
    }

	// Here
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 1; j < NY - 1; j++){
            U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j+1,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k+1,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }
    }

	// There
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 1; j < NY - 1; j++){
            U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j+1,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }
    }

    // Bottom Here Corner
    j = 0;
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j+1,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k+1,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
    }

    // Bottom There Corner
    j = 0;
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j+1,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
    }

    // Top Here Corner
    j = NY - 1;
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Top[TOP(i,j,k)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k+1,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
    }

    // Top There Corner
    j = NY - 1;
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Top[TOP(i,j,k)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U.Pres[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U.Pres[LU(i,j,k,0)] - U.Pres[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
    }

    if (Rango == 0){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){	
				U.Diffusive[LU(0,j,k,0)] = 0.0;
			}
		}
    }
    else if (Rango == Procesos - 1){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){	
				U.Diffusive[LU(NX,j,k,0)] = 0.0;
			}
		}
    }

}

// Function to calculate the diffusion term of V Velocity
void CFD_Solver::Get_DiffusionV(Mesher MESH){
int i, j, k;

    // Center
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i+1,j,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k+1,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
		}
	}

    // Here
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i+1,j,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k+1,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}
	}

    // There
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i+1,j,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}
	}

    if (Rango == 0){

        i = 0;
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i+1,j,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k+1,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
		}

        // Here
        k = 0;
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i+1,j,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k+1,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

        // There
        k = NZ - 1;
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i+1,j,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

    }
    else if (Rango == Procesos - 1){

        i = NX - 1;
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Right[RIGHT(i,j,k)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k+1,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
		}

        // Here
        k = 0;
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Right[RIGHT(i,j,k)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k+1,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

        // There
        k = NZ - 1;
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Right[RIGHT(i,j,k)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V.Pres[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V.Pres[LV(i,j,k,0)] - V.Pres[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

	}

}

// Function to calculate the diffusion term of W Velocity
void CFD_Solver::Get_DiffusionW(Mesher MESH){
int i, j, k;

    // Center
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i+1,j,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j+1,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
			}
		}
    }

    // Bottom
    j = 0;
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 1; k < NZ; k++){
			W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i+1,j,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j+1,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}
    }

    // Top
    j = NY - 1;
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 1; k < NZ; k++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i+1,j,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Top[TOP(i,j,k)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}
    }

    if (Rango == 0){

        i = 0;
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i+1,j,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j+1,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
			}
		}

        // Bottom
        j = 0;
		for(k = 1; k < NZ; k++){
			W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i+1,j,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j+1,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}
    
        // Top
        j = NY - 1;
		for(k = 1; k < NZ; k++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i+1,j,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Top[TOP(i,j,k)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}

    }
    else if(Rango == Procesos - 1){

        i = NX - 1;
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Right[RIGHT(i,j,k)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j+1,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
			}
		}

        // Bottom
        j = 0;
		for(k = 1; k < NZ; k++){
			W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Right[RIGHT(i,j,k)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j+1,k,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}

        // Top
        j = NY - 1;
		for(k = 1; k < NZ; k++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Right[RIGHT(i,j,k)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Top[TOP(i,j,k)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,0)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] - W.Pres[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W.Pres[LW(i,j,k,0)] - W.Pres[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}

    }

}
