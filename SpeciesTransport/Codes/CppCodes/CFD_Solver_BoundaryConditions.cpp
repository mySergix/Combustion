//------------------------------------------------------------------------------------------------//
//                            CPP FILE FOR BOUNDARY CONDITIONS SETUP                              //
//------------------------------------------------------------------------------------------------//

// Function to get the boundary conditions of each property
void CFD_Solver::Get_BoundaryConditions(){
int i, j, k;

    // Parte Top y Bottom
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (k = - Halo; k < NZ + Halo; k++){

            // Parte Bottom
            Density.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // Density
            U.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // U Velocity
            V.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // V Velocity
            W.Bottom[TOP_BOT(i,0,k,0)] = 0.0; // V Velocity

            //Parte Top
            Density.Top[TOP_BOT(i,NY,k,0)] = 0.0; //Density
            U.Top[TOP_BOT(i,NY,k,0)] = 0.0; // U Velocity
            V.Top[TOP_BOT(i,NY,k,0)] = 0.0; // V Velocity
            W.Top[TOP_BOT(i,NY,k,0)] = 0.0; // W Velocity

        }
    }

    // Parte Here y There
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){

            // Parte Here
            Density.Here[HER_THE(i,j,0,0)] = 0.0; // Density
            

            // Parte There 
            Density.There[HER_THE(i,j,NZ,0)] = 0.0; // Density
            U.There[HER_THE(i,j,NZ,0)] = 0.0; // Velocity U
            V.There[HER_THE(i,j,NZ,0)] = 0.0; // Velocity V
            W.There[HER_THE(i,j,NZ,0)] = 0.0; // Velocity W

        }
    }

    // Parte Left
    if (Rango == 0){

        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){

                Density.Left[LEF_RIG(0,j,k,0)] = 1.0; // Density
                U.Left[LEF_RIG(0,j,k,0)] = 10.0; // Velocity U
                V.Left[LEF_RIG(0,j,k,0)] = 0.0; // Velocity V
                W.Left[LEF_RIG(0,j,k,0)] = 0.0; // Velocity W

            }
        }

    }

}

// Function to calculate and update periodic conditions
void CFD_Solver::Get_PeriodicConditions(Mesher MESH, Property_Struct &PropertyName){
int i, j, k;

    // Here
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < 0; k++){
                PropertyName.Pres[LM(i,j,k,0)] = PropertyName.Pres[LM(i,j,NZ+k,0)];
            }
        }
    }

    // There
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = NZ; k < NZ + Halo; k++){
                PropertyName.Pres[LM(i,j,k,0)] = PropertyName.Pres[LM(i,j,k-NZ,0)];
            }
        }
    } 

    // Walls
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            PropertyName.Here[HER_THE(i,j,0,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,-1,2)] + MESH.Node_Mesh[LM(i,j,0,2)]), 0.50 * (W.Pres[LM(i,j,-1,0)] + W.Pres[LM(i,j,0,0)]), MESH.Node_Mesh[LM(i,j,-2,2)], PropertyName.Pres[LM(i,j,-2,0)], MESH.Node_Mesh[LM(i,j,-1,2)], PropertyName.Pres[LM(i,j,-1,0)], MESH.Node_Mesh[LM(i,j,0,2)], PropertyName.Pres[LM(i,j,0,0)], MESH.Node_Mesh[LM(i,j,1,2)], PropertyName.Pres[LM(i,j,1,0)]);
            PropertyName.There[HER_THE(i,j,NZ,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,NZ-1,2)] + MESH.Node_Mesh[LM(i,j,NZ,2)]), 0.50 * (W.Pres[LM(i,j,NZ-1,0)] + W.Pres[LM(i,j,NZ,0)]), MESH.Node_Mesh[LM(i,j,NZ-2,2)], PropertyName.Pres[LM(i,j,NZ-2,0)], MESH.Node_Mesh[LM(i,j,NZ-1,2)], PropertyName.Pres[LM(i,j,NZ-1,0)], MESH.Node_Mesh[LM(i,j,NZ,2)], PropertyName.Pres[LM(i,j,NZ,0)], MESH.Node_Mesh[LM(i,j,NZ+1,2)], PropertyName.Pres[LM(i,j,NZ+1,0)]);
        }
    }

}

// Function to update continuously the boundary conditions
void CFD_Solver::Update_BoundaryConditions(Mesher MESH){
int i, j, k;

    Get_PeriodicConditions(MESH, U); // Periodic Conditions Velocity U
    Get_PeriodicConditions(MESH, V); // Periodic Conditions Velocity V
    Get_PeriodicConditions(MESH, W); // Periodic Conditions Velocity W

    if (Rango == Procesos - 1){

        // Right Side
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Density.Right[LEF_RIG(0,j,k,0)] = Density.Pres[LM(NX-1,j,k,0)]; // Density
                U.Right[LEF_RIG(0,j,k,0)] = U.Pres[LM(NX-1,j,k,0)]; // Velocity U
                V.Right[LEF_RIG(0,j,k,0)] = V.Pres[LM(NX-1,j,k,0)]; // Velocity V
                W.Right[LEF_RIG(0,j,k,0)] = W.Pres[LM(NX-1,j,k,0)]; // Velocity W
            }
        }
        
        // Pressure Right Side
        for (i = NX; i < NX + Halo; i++){
            for (j = - Halo; j < NY + Halo; j++){
                for (k = - Halo; k < NZ + Halo; k++){
                    Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(NX-1,j,k,0)]; // Pressure
                }
            }
        }

    }

    if (Rango == 0){

        // Pressure Left Side
        for (i = - Halo; i < 0; i++){
            for (j = - Halo; j < NY + Halo; j++){
                for (k = - Halo; k < NZ + Halo; k++){
                    Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(0,j,k,0)]; // Pressure
                }
            }
        }

    }

    // Bottom 
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < 0; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(i,0,k,0)];
            }
        }
    }

    // Pressure Top
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[LM(i,NY-1,k,0)];
            }
        }
    }

    // Pressure Here
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < 0; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[i,j,NX + k,0];
            }
        }
    }

    // Pressure There
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = NX; k < NX + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] = Pressure.Pres[i,j,k - NX,0];
            }
        }
    }

}

// Function for the application of the boundary condition to the walls of the volumes
void CFD_Solver::ApplyBoundaries(Property_Struct &PropertyName){
int i, j, k;

    // Left Side
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){               
                PropertyName.Wall_U[LMU(0,j,k,0)] = PropertyName.Left[LEF_RIG(0,j,k,0)];
            }
        }
    }

    // Right Side
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){               
                PropertyName.Wall_U[LMU(NX,j,k,0)] = PropertyName.Right[LEF_RIG(NX,j,k,0)];
            }
        }
    }

    // Top and Bottom Sides
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            // Bottom Part
            PropertyName.Wall_V[LMV(i,0,k,0)] = PropertyName.Bottom[TOP_BOT(i,0,k,0)];

            // Top Part
            PropertyName.Wall_V[LMV(i,NY,k,0)] = PropertyName.Top[TOP_BOT(i,NY,k,0)];
        }
    }

    // Here and There Sides
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            // Here Part
            PropertyName.Wall_W[LMW(i,j,0,0)] = PropertyName.Here[HER_THE(i,j,0,0)];

            // There Part
            PropertyName.Wall_W[LMW(i,j,NZ,0)] = PropertyName.There[HER_THE(i,j,NZ,0)];
        }
    }

}