//------------------------------------------------------------------------------------------------//
//                            CPP FILE FOR BOUNDARY CONDITIONS SETUP                              //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Boundary Conditions of the simulation
void CFD_Solver::Get_InitialBoundaryConditions(Mesher MESH){
int i, j, k;

    // Velocities
    double m = 1.7 + 0.50 * pow(gamma_geometry, -1.4);
    double n = 2 + 0.3 * (gamma_geometry - 1/3);
    if (Rango == 0){ // Left
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Left[LEFT(0,j,k)] = 0.0;//w_av * ((m+1)/m) * ((n+1)/n) * (1.0 - pow(MESH.MP[LP(0,j,k,1)]/Ydominio, n)) * (1.0 - pow(MESH.MP[LP(0,j,k,2)]/Zdominio, m));
                V.Left[LEFT(0,j,k)] = 0.0;
                W.Left[LEFT(0,j,k)] = 0.0;
            }
        }
    }
     
    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            U.Bottom[BOTTOM(i,0,k)] = 0.0;
            V.Bottom[BOTTOM(i,0,k)] = 0.0;
            W.Bottom[BOTTOM(i,0,k)] = 0.0;
        }
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            U.Top[TOP(i,NY,k)] = 0.0;
            V.Top[TOP(i,NY,k)] = 0.0;
            W.Top[TOP(i,NY,k)] = 0.0;
        }
    }
    
     // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.Here[HERE(i,j,0)] = 0.0;
            V.Here[HERE(i,j,0)] = 0.0;
            W.Here[HERE(i,j,0)] = 0.0;
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.There[THERE(i,j,0)] = 0.0;
            V.There[THERE(i,j,0)] = 0.0;
            W.There[THERE(i,j,0)] = 0.0;
        }
    }


    // Temperaturas
    if (Rango == 0){ // Left
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                T.Left[LEFT(0,j,k)] = To;
            }
        }
    }

    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            T.Bottom[BOTTOM(i,0,k)] = T_water;
        }
    }

}

// Function to calculate the periodic boundary conditions
void CFD_Solver::Get_UpdateBoundaryConditions(Mesher MESH){
int i, j, k;

    // Velocities
    double m = 1.7 + 0.50 * pow(gamma_geometry, -1.4);
    double n = 2 + 0.3 * (gamma_geometry - 1/3);
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Right[RIGHT(NX,j,k)] = U.Pres[LU(NX-1,j,k,0)];
                V.Right[RIGHT(NX,j,k)] = V.Pres[LV(NX-1,j,k,0)];
                W.Right[RIGHT(NX,j,k)] = W.Pres[LW(NX-1,j,k,0)];
            }
        }
    }  

   
    // Temperatures

    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                T.Right[RIGHT(NX,j,k)] = T.Pres[LP(NX-1,j,k,0)];
            }
        }
    }  

    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            T.Here[HERE(i,j,0)] = T.Pres[LP(i,j,0,0)];
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            T.There[THERE(i,j,0)] = T.Pres[LP(i,j,NZ-1,0)];
        }
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            T.Top[TOP(i,0,k)] = T.Pres[LP(i,NY-1,k,0)];
        }
    }

}

// Function to set all the corresponding velocity halos
void CFD_Solver::Get_StaticHalos(){
int i, j, k;

    // Velocities
    if (Rango == 0){ // Left
        for (i = - Halo; i < Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    U.Pres[LU(i,j,k,0)] = U.Left[LEFT(0,j,k)];
                    V.Pres[LV(i,j,k,0)] = V.Left[LEFT(0,j,k)];
                    W.Pres[LW(i,j,k,0)] = W.Left[LEFT(0,j,k)];
                }
            }
        }  
    }
    

    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = - Halo; j < 0; j++){
            for (k = 0; k < NZ; k++){
            U.Pres[LU(i,j,k,0)] = U.Bottom[BOTTOM(i,0,k)];
            V.Pres[LV(i,j,k,0)] = V.Bottom[BOTTOM(i,0,k)];
            W.Pres[LV(i,j,k,0)] = W.Bottom[BOTTOM(i,0,k)];
            }
        }   
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = 0; k < NZ; k++){
            U.Pres[LU(i,j,k,0)] = U.Top[TOP(i,NY,k)];
            V.Pres[LV(i,j,k,0)] = V.Top[TOP(i,NY,k)];
            W.Pres[LW(i,j,k,0)] = W.Top[TOP(i,NY,k)];
            }
        }  
    }

    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k < 0; k++){
                U.Pres[LU(i,j,k,0)] = U.Here[HERE(i,j,0)];
                V.Pres[LV(i,j,k,0)] = V.Here[HERE(i,j,0)];
                W.Pres[LW(i,j,k,0)] = W.Here[HERE(i,j,0)];
            }  
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo; k++){
                U.Pres[LU(i,j,k,0)] = U.There[THERE(i,j,0)];
                V.Pres[LV(i,j,k,0)] = V.There[THERE(i,j,0)];
                W.Pres[LW(i,j,k,0)] = W.There[THERE(i,j,0)];
            }
            
        }
    }


    // Temperatures
    if (Rango == 0){ // Left
        for (i = - Halo; i < Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    T.Pres[LP(i,j,k,0)] = T.Left[LEFT(0,j,k)];   
                }
            }
        }  
    }
     
}

// Function to update the boundary conditions
void CFD_Solver::Get_UpdateHalos(){
int i, j, k;

    // Velocities

    if (Rango == Procesos - 1){
        for (i = Fx[Rango]; i < Fx[Rango] + Halo; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    U.Pres[LU(i,j,k,0)] = U.Right[RIGHT(NX,j,k)];
                    V.Pres[LV(i,j,k,0)] = V.Right[RIGHT(NX,j,k)];
                    W.Pres[LW(i,j,k,0)] = W.Right[RIGHT(NX,j,k)];
                }
            }
        }    
    }   

    
    // Temperatures

    if (Rango == Procesos - 1){
        for (i = Fx[Rango]; i < Fx[Rango] + Halo; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    T.Pres[LP(i,j,k,0)] = T.Right[RIGHT(NX,j,k)];
                }
            }
        }    
    }  

    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = - Halo; j < 0; j++){
            for (k = 0; k < NZ; k++){
                T.Pres[LP(i,j,k,0)] = T.Bottom[BOTTOM(i,0,k)];
            }
        }   
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = 0; k < NZ; k++){
                T.Pres[LP(i,j,k,0)] = T.Top[TOP(i,NY,k)];
            }
        }  
    }

    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k < 0; k++){
                T.Pres[LP(i,j,k,0)] = T.Here[HERE(i,j,0)];
            }  
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo; k++){
                T.Pres[LP(i,j,k,0)] = T.There[THERE(i,j,0)];
            }
        }
    }

}