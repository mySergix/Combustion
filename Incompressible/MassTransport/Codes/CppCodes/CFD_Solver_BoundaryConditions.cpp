//------------------------------------------------------------------------------------------------//
//                            CPP FILE FOR BOUNDARY CONDITIONS SETUP                              //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Boundary Conditions of the simulation
void CFD_Solver::Get_InitialBoundaryConditions(){
int i, j, k;

    if (Rango == 0){ // Left
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Left[LEFT(0,j,k)] = 0.0;
                V.Left[LEFT(0,j,k)] = 0.0;
                W.Left[LEFT(0,j,k)] = 0.0;
            }
        }
    }
    else if (Rango == Procesos - 1){ // Right
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Right[RIGHT(NX,j,k)] = 0.0;
                V.Right[RIGHT(NX,j,k)] = 0.0;
                W.Right[RIGHT(NX,j,k)] = 0.0;
            }
        }
    }

    // Bottom and Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            U.Bottom[BOTTOM(i,0,k)] = 0.0;
            V.Bottom[BOTTOM(i,0,k)] = 0.0;
            W.Bottom[BOTTOM(i,0,k)] = 0.0;

            U.Top[TOP(i,NY,k)] = Uref;
            V.Top[TOP(i,NY,k)] = 0.0;
            W.Top[TOP(i,NY,k)] = 0.0;
        }
    }

}

// Function to calculate the periodic boundary conditions
void CFD_Solver::Get_PeriodicBoundaryConditions(){
int i, j, k;

    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.Here[HERE(i,j,0)] = 0.50 * (U.Pres[LU(i,j,NZ-1,0)] + U.Pres[LU(i+1,j,NZ-1,0)]);
            V.Here[HERE(i,j,0)] = 0.50 * (V.Pres[LV(i,j,NZ-1,0)] + V.Pres[LV(i,j+1,NZ-1,0)]);
            W.Here[HERE(i,j,0)] = W.Pres[LW(i,j,NZ-1,0)];
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.There[THERE(i,j,0)] = 0.50 * (U.Pres[LU(i,j,0,0)] + U.Pres[LU(i+1,j,0,0)]);
            V.There[THERE(i,j,0)] = 0.50 * (V.Pres[LV(i,j,0,0)] + V.Pres[LV(i,j+1,0,0)]);
            W.There[THERE(i,j,0)] = W.Pres[LW(i,j,1,0)];
        }
    }

}

// Function to set all the corresponding velocity halos
void CFD_Solver::Get_StaticHalos(){
int i, j, k;

    // Velocity U
    if (Rango == 0){
        for (i = - Halo; i < Ix[Rango] + 1; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    U.Pres[LU(i,j,k,0)] = U.Left[LEFT(i,j,k)];
                }
            }
        }
    }
    else if (Rango == Procesos - 1){
        for (i = Fx[Rango]; i < Fx[Rango] + Halo + 1; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    U.Pres[LU(i,j,k,0)] = U.Right[RIGHT(i,j,k)];
                }
            }
        }
    }


}