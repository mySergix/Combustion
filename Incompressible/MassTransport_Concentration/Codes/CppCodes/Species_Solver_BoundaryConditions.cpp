//------------------------------------------------------------------------------------------------//
//                        CPP FILE FOR SPECIES BOUNDARY CONDITIONS SETUP                          //
//------------------------------------------------------------------------------------------------//

// Function to set the initial fields
void Species_Solver::Get_InitialConditionsSpecies(){
int i, j, k;

	// Temperature T
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){	
				Species[0].C_Pres[LP(i,j,k,0)] = 0.0;
				Species[0].C_Fut[LP(i,j,k,0)] = 0.0;
			}
		}
	}
	
}

// Function to set the initial boundary conditions of the species
void Species_Solver::Get_InitialBoundaryConditionsSpecies(){
int i, j, k;

    // Here (Inlet)
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            Species[0].Here[HERE(i,j,0)] = Co;
        }
    }

}

// Function to update the boundary conditions of the species
void Species_Solver::Get_UpdateBoundaryConditionsSpecies(Mesher MESH, double *Tmatrix){
int i, j, k;

    // Bottom (Water)
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            if (MESH.MP[LP(i,0,k,0)] >= 0.09 && MESH.MP[LP(i,0,k,0)] <= Xdominio - 0.09){
                Species[0].Bottom[BOTTOM(i,0,k)] = 1.0 + ((T_water - Tmatrix[LP(i,0,k,0)]) / (hfg * D_AB));
            }
            else{
                Species[0].Bottom[BOTTOM(i,0,k)] = Species[0].C_Pres[LP(i,0,k,0)];
            }       
        }
    }

    // There (Outlet) (Null Gradients)
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            Species[0].There[THERE(i,j,0)] = Species[0].C_Pres[LP(i,j,NZ-1,0)];
        }
    }

    if (Rango == 0){ // Left
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[0].Left[LEFT(0,j,k)] = Species[0].C_Pres[LP(0,j,k,0)];
            }
        }
    }
    else if (Rango == Procesos - 1){ // Right
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[0].Right[RIGHT(NX,j,k)] = Species[0].C_Pres[LP(NX-1,j,k,0)];
            }
        }
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            Species[0].Top[TOP(i,NY,k)] = Species[0].C_Pres[LP(i,NY-1,k,0)];
        }
    }

}

// Function to get all the initial halos of the species boundary conditions
void Species_Solver::Get_InitialHalosSpecies(){
int i, j, k;

    // Here (Inlet)
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = - HP; k < 0; k++){
                Species[0].C_Pres[LP(i,j,k,0)] = Species[0].Here[HERE(i,j,k)];
            }
        }
    }

}

// Function to update all the halos of the species boundary conditions
void Species_Solver::Get_UpdateHalosSpecies(){
int i, j, k;

    // Bottom (Water)
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = - HP; j < 0; j++){
            for (k = 0; k < NZ; k++){
                Species[0].C_Pres[LP(i,j,k,0)] = Species[0].Bottom[BOTTOM(i,0,k)];      
            }
        }
    }

    // There (Outlet) (Null Gradients)
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + HP; k++){
                Species[0].C_Pres[LP(i,j,k,0)] = Species[0].There[THERE(i,j,0)];
            } 
        }
    }

    if (Rango == 0){ // Left
        for (i = - HP; i < Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    Species[0].C_Pres[LP(i,j,k,0)] = Species[0].Left[LEFT(0,j,k)];
                }
            }
        }  
    }
    else if (Rango == Procesos - 1){ // Right
        for (i = NX; i < NX + HP; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    Species[0].C_Pres[LP(i,j,k,0)] = Species[0].Right[RIGHT(0,j,k)];
                }
            }
        }
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + HP; j++){
            for (k = 0; k < NZ; k++){
                Species[0].C_Pres[LP(i,j,k,0)] = Species[0].Top[TOP(i,NY,k)];
            }
        }
        
    }

}