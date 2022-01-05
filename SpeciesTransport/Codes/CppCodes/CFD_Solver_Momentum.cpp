//------------------------------------------------------------------------------------------------//
//                        CPP FILE FOR N-S MOMENTUM EQUATION CALCULATIONS                         //
//------------------------------------------------------------------------------------------------//

// Calculation of Velocities at the Walls of the Volumes
void CFD_Solver::Get_WallsVelocities(Mesher MESH){
int i, j, k;

    // MU Walls
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], U.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], U.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], U.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], U.Pres[LM(i+1,j,k,0)]);
                V.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], V.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], V.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], V.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], V.Pres[LM(i+1,j,k,0)]);
                W.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], W.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], W.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], W.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], W.Pres[LM(i+1,j,k,0)]);
            }
        }
    }

    ApplyBoundaries(U);  

    // MV Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], U.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], U.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], U.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], U.Pres[LM(i,j+1,k,0)]);
                V.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], V.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], V.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], V.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], V.Pres[LM(i,j+1,k,0)]);
                W.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], W.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], W.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], W.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], W.Pres[LM(i,j+1,k,0)]);
            }
        }
    }

    ApplyBoundaries(V);  

    // MW Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 1; k < NZ; k++){
                U.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], U.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], U.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], U.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], U.Pres[LM(i,j,k+1,0)]);
                V.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], V.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], V.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], V.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], V.Pres[LM(i,j,k+1,0)]);
                W.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], W.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], W.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], W.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], W.Pres[LM(i,j,k+1,0)]);
            }
        }
    }

    ApplyBoundaries(W);  

}

// Function to compute the Convective Term of a Property
void CFD_Solver::Get_ConvectiveTerm(Mesher MESH, Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Convective[LM(i,j,k,0)] = (1.0/MESH.Vol[LM(i,j,k,0)]) * (
                                                     + MESH.Surf[LM(i,j,k,0)] * Density.Wall_U[LMU(i+1,j,k,0)] * U.Wall_U[LMU(i+1,j,k,0)] * PropertyName.Wall_U[LMU(i+1,j,k,0)]
                                                     - MESH.Surf[LM(i,j,k,0)] * Density.Wall_U[LMU(i,j,k,0)] * U.Wall_U[LMU(i,j,k,0)] * PropertyName.Wall_U[LMU(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,1)] * Density.Wall_V[LMV(i,j+1,k,0)] * V.Wall_V[LMV(i,j+1,k,0)] * PropertyName.Wall_V[LMV(i,j+1,k,0)]
                                                     - MESH.Surf[LM(i,j,k,1)] * Density.Wall_V[LMV(i,j,k,0)] * V.Wall_V[LMV(i,j,k,0)] * PropertyName.Wall_V[LMV(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,2)] * Density.Wall_W[LMW(i,j,k+1,0)] * W.Wall_W[LMW(i,j,k+1,0)] * PropertyName.Wall_W[LMW(i,j,k+1,0)]
                                                     - MESH.Surf[LM(i,j,k,2)] * Density.Wall_W[LMW(i,j,k,0)] * W.Wall_W[LMW(i,j,k,0)] * PropertyName.Wall_W[LMW(i,j,k,0)]
                                                     );
            }
        }
    }

}

//Function to calculate the pressure
void CFD_Solver::Get_Pressure(){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Pressure.Pres[LM(i,j,k,0)] = Density.Pres[LM(i,j,k,0)] * 287.0 * 298.15;
            }
        }
    }

}

// Function to calculate the pressure gradient
void CFD_Solver::Get_PressureGradient(Mesher MESH){
int i, j, k;
double Px_1, Px_2, Py_1, Py_2, Pz_1, Pz_2;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Px_1 = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], Pressure.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], Pressure.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], Pressure.Pres[LM(i+1,j,k,0)]);
                Px_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,0)] + MESH.Node_Mesh[LM(i + 1,j,k,0)]), 0.50 * (U.Pres[LM(i,j,k,0)] + U.Pres[LM(i + 1,j,k,0)]), MESH.Node_Mesh[LM(i-1,j,k,0)], Pressure.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], Pressure.Pres[LM(i+1,j,k,0)], MESH.Node_Mesh[LM(i+2,j,k,0)], Pressure.Pres[LM(i+2,j,k,0)]);
                
                Py_1 = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], Pressure.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], Pressure.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], Pressure.Pres[LM(i,j+1,k,0)]);
                Py_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,1)] + MESH.Node_Mesh[LM(i,j+1,k,1)]), 0.50 * (V.Pres[LM(i,j,k,0)] + V.Pres[LM(i,j+1,k,0)]), MESH.Node_Mesh[LM(i,j-1,k,1)], Pressure.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], Pressure.Pres[LM(i,j+1,k,0)], MESH.Node_Mesh[LM(i,j+2,k,1)], Pressure.Pres[LM(i,j+2,k,0)]);
                
                Pz_1 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], Pressure.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], Pressure.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], Pressure.Pres[LM(i,j,k+1,0)]);
                Pz_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,2)] + MESH.Node_Mesh[LM(i,j,k+1,2)]), 0.50 * (W.Pres[LM(i,j,k,0)] + W.Pres[LM(i,j,k+1,0)]), MESH.Node_Mesh[LM(i,j,k-1,2)], Pressure.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], Pressure.Pres[LM(i,j,k+1,0)], MESH.Node_Mesh[LM(i,j,k+2,2)], Pressure.Pres[LM(i,j,k+2,0)]);

                Pressure.Gradient_X[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,0)]) * (Px_2 - Px_1);
                Pressure.Gradient_Y[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,1)]) * (Py_2 - Py_1);
                Pressure.Gradient_Z[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,2)]) * (Pz_2 - Pz_1);
            }
        }
    }

}

// Function to calculate the divergence of the field
void CFD_Solver::Get_Divergence(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Divergence[LM(i,j,k,0)] = (1.0 / MESH.Vol[LM(i,j,k,0)]) * (
                                        - MESH.Surf[LM(i,j,k,0)] * U.Wall_U[LMU(i,j,k,0)]
                                        + MESH.Surf[LM(i,j,k,0)] * U.Wall_U[LMU(i+1,j,k,0)]
                                        - MESH.Surf[LM(i,j,k,1)] * V.Wall_V[LMV(i,j,k,0)]
                                        + MESH.Surf[LM(i,j,k,1)] * V.Wall_V[LMV(i,j+1,k,0)]
                                        - MESH.Surf[LM(i,j,k,2)] * W.Wall_W[LMW(i,j,k,0)]
                                        + MESH.Surf[LM(i,j,k,2)] * W.Wall_W[LMW(i,j,k+1,0)]
                                        );
            }
        }
    }

}

// Function to calculate the velociy gradients at the walls
void CFD_Solver::Get_VelocityGradients(Mesher MESH, Property &PropertyName){
int i, j, k;

    // Gradients in X direction

        // Walls U
        for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    PropertyName.Wall_U_Diff_X[LMU(i,j,k,0)] = (PropertyName.Pres[LM(i,j,k,0)] - PropertyName.Pres[LM(i-1,j,k,0)]) / (0.50 * (MESH.DeltaP[LM(i,j,k,0)] + MESH.DeltaP[LM(i-1,j,k,0)]));
                }
            }
        }

        // Walls V
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 0; j < NY + 1; j++){
                for (k = 0; k < NZ; k++){
                    PropertyName.Wall_V_Diff_X[LMV(i,j,k,0)] = 0.50 * (PropertyName.Wall_V[LMV(i+1,j,k,0)] - PropertyName.Wall_V[LMV(i-1,j,k,0)]) / (MESH.DeltaP[LM(i,j,k,0)]);
                }
            }
        }

        // Walls W
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    PropertyName.Wall_W_Diff_X[LMW(i,j,k,0)] = 0.50 * (PropertyName.Wall_W[LMW(i+1,j,k,0)] - PropertyName.Wall_W[LMW(i-1,j,k,0)]) / (MESH.DeltaP[LM(i,j,k,0)]);
                }
            }
        }

    // Gradients in Y direction

        // Walls U
        for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for (j = 1; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    PropertyName.Wall_U_Diff_Y[LMU(i,j,k,0)] = 0.50 * (PropertyName.Wall_U[LMU(i,j+1,k,0)] - PropertyName.Wall_U[LMU(i,j-1,k,0)]) / (MESH.DeltaP[LM(i,j,k,1)]);
                    }
            }
        }

        // Walls V
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 1; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    // Parte Top
                    PropertyName.Wall_V_Diff_Y[LMV(i,0,k,0)] = 0.0;

                    // Parte Bottom
                    PropertyName.Wall_V_Diff_Y[LMV(i,NY,k,0)] = 0.0;

                    // Resto
                    PropertyName.Wall_V_Diff_Y[LMV(i,j,k,0)] = (PropertyName.Pres[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j-1,k,0)]) / (0.50 * (MESH.DeltaP[LM(i,j,k,1)] + MESH.DeltaP[LM(i,j-1,k,1)]));
                }
            }
        }

        // Walls W
        for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for (j = 1; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    PropertyName.Wall_W_Diff_Y[LMW(i,j,k,0)] = 0.50 * (PropertyName.Wall_W[LMW(i,j+1,k,0)] - PropertyName.Wall_W[LMW(i,j-1,k,0)]) / (MESH.DeltaP[LM(i,j,k,1)]);
                    }
            }
        }

    // Gradients in Z direction

        // Walls U
        for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    PropertyName.Wall_U_Diff_Z[LMU(i,j,k,0)] = 0.50 * (PropertyName.Wall_U[LMU(i,j,k+1,0)] - PropertyName.Wall_U[LMU(i,j,k-1,0)]) / (MESH.DeltaP[LM(i,j,k,2)]);
                }
            }
        }

        // Walls V
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ + 1; k++){
                    PropertyName.Wall_V_Diff_Z[LMV(i,j,k,0)] = 0.50 * (PropertyName.Wall_V[LMV(i,j,k+1,0)] - PropertyName.Wall_V[LMV(i,j,k-1,0)]) / (MESH.DeltaP[LM(i,j,k,2)]);
                }
            }
        }

        // Walls W
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    PropertyName.Wall_W_Diff_Z[LMW(i,j,k,0)] = (PropertyName.Pres[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j,k-1,0)]) / (0.50 * (MESH.DeltaP[LM(i,j,k,2)] + MESH.DeltaP[LM(i,j,k-1,2)]));
                }
            }
        }
}

// Function to calculate the viscous stresses of the fluid
void CFD_Solver::Get_ViscousStresses(Mesher MESH){
int i, j, k;

    // Tau XX Walls U
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_xx_X[LMU(i,j,k,0)] =  - (2.0/3.0) * mu * 0.50 * (Stress.Divergence[LM(i-1,j,k,0)] + Stress.Divergence[LM(i,j,k,0)]) + 2.0 * mu * U.Wall_U_Diff_X[LMU(i,j,k,0)];
            }
        }
    }

    // Tau YY Walls V
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY + 1; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_yy_Y[LMV(i,j,k,0)] =  - (2.0/3.0) * mu * 0.50 * (Stress.Divergence[LM(i,j-1,k,0)] + Stress.Divergence[LM(i,j,k,0)]) + 2.0 * mu * V.Wall_V_Diff_Y[LMV(i,j,k,0)];
            }
        }
    }

    // Tau ZZ Walls W
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                Stress.Tau_zz_Z[LMW(i,j,k,0)] =  - (2.0/3.0) * mu * 0.50 * (Stress.Divergence[LM(i,j,k-1,0)] + Stress.Divergence[LM(i,j,k,0)]) + 2.0 * mu * W.Wall_W_Diff_Z[LMW(i,j,k,0)];
            }
        }
    }

    // Tau XY Walls V
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY + 1; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_yx_Y[LMV(i,j,k,0)] =  mu * (U.Wall_V_Diff_Y[LMV(i,j,k,0)] + V.Wall_V_Diff_X[LMV(i,j,k,0)]);
            }
        }
    }

    // Tau ZX Walls W
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                Stress.Tau_zx_Z[LMW(i,j,k,0)] =  mu * (W.Wall_W_Diff_X[LMW(i,j,k,0)] + U.Wall_W_Diff_Z[LMW(i,j,k,0)]);
            }
        }
    }

    // Tau XY Walls U
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_xy_X[LMU(i,j,k,0)] =  mu * (U.Wall_U_Diff_Y[LMU(i,j,k,0)] + V.Wall_U_Diff_X[LMU(i,j,k,0)]);
            }
        }
    }

    // Tau ZY Walls W
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                Stress.Tau_zy_Z[LMW(i,j,k,0)] =  mu * (V.Wall_W_Diff_Z[LMW(i,j,k,0)] + W.Wall_W_Diff_Y[LMW(i,j,k,0)]);
            }
        }
    }

    // Tau YZ Walls V
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY + 1; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_yz_Y[LMV(i,j,k,0)] =  mu * (V.Wall_V_Diff_Z[LMV(i,j,k,0)] + W.Wall_V_Diff_Y[LMV(i,j,k,0)]);
            }
        }
    }

    // Tau XZ Walls U
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_xz_X[LMU(i,j,k,0)] =  mu * (W.Wall_U_Diff_X[LMU(i,j,k,0)] + U.Wall_U_Diff_Z[LMU(i,j,k,0)]);
            }
        }
    }

}

// Function to calculate the diffusive term of each velocity
void CFD_Solver::Get_DiffusiveTermVelocity(Mesher MESH, Property &PropertyName, double *Wall_U_Term, double *Wall_V_Term, double *Wall_W_Term){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Diffusive[LM(i,j,k,0)] = (1.0/MESH.Vol[LM(i,j,k,0)]) * (
                                                     - MESH.Surf[LM(i,j,k,0)] * Wall_U_Term[LMU(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,0)] * Wall_U_Term[LMU(i+1,j,k,0)]
                                                     - MESH.Surf[LM(i,j,k,1)] * Wall_V_Term[LMV(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,1)] * Wall_V_Term[LMV(i,j+1,k,0)]
                                                     - MESH.Surf[LM(i,j,k,2)] * Wall_W_Term[LMW(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,2)] * Wall_W_Term[LMW(i,j,k+1,0)]
                                                     );
            }
        }
    }

}

// Function to calculate the step contribution to each velocity
void CFD_Solver::Get_StepContribution_Velocity(Property &PropertyName, double *PressureGradient){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.ContributionPres[LM(i,j,k,0)] = (1.0 / Density.Pres[LM(i,j,k,0)]) * (- PropertyName.Convective[LM(i,j,k,0)] - PressureGradient[LM(i,j,k,0)] + PropertyName.Diffusive[LM(i,j,k,0)]) + PropertyName.Gravity;
            }
        }
    }

}