//------------------------------------------------------------------------------------------------//
//                        CPP FILE FOR N-S MOMENTUM EQUATION CALCULATIONS                         //
//------------------------------------------------------------------------------------------------//

// Function to compute the Convective Term of a Property
void CFD_Solver::Get_ConvectiveTerm_Property(Mesher MESH, Property_Struct &PropertyName){
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
                Pressure.Pres[LM(i,j,k,0)] = Density.Pres[LM(i,j,k,0)] * R_ideal * T_Pres[LM(i,j,k,0)];
            }
        }
    }

}

// Function to calculate the pressure gradient
void CFD_Solver::Get_PressureGradient(Mesher MESH){
int i, j, k;
double P_1, P_2;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){

                P_1 = CS(0.50 * (MESH.Node_Mesh[LM(i - 1,j,k,0)] + MESH.Node_Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i-2,j,k,0)], Pressure.Pres[LM(i-2,j,k,0)], MESH.Node_Mesh[LM(i-1,j,k,0)], Pressure.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], Pressure.Pres[LM(i+1,j,k,0)]);
                P_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,0)] + MESH.Node_Mesh[LM(i + 1,j,k,0)]), 0.50 * (U.Pres[LM(i,j,k,0)] + U.Pres[LM(i + 1,j,k,0)]), MESH.Node_Mesh[LM(i-1,j,k,0)], Pressure.Pres[LM(i-1,j,k,0)], MESH.Node_Mesh[LM(i,j,k,0)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i+1,j,k,0)], Pressure.Pres[LM(i+1,j,k,0)], MESH.Node_Mesh[LM(i+2,j,k,0)], Pressure.Pres[LM(i+2,j,k,0)]);
                
                Pressure.Gradient_X[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,0)]) * (P_2 - P_1);

                P_1 = CS(0.50 * (MESH.Node_Mesh[LM(i,j-1,k,1)] + MESH.Node_Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j-2,k,1)], Pressure.Pres[LM(i,j-2,k,0)], MESH.Node_Mesh[LM(i,j-1,k,1)], Pressure.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], Pressure.Pres[LM(i,j+1,k,0)]);
                P_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,1)] + MESH.Node_Mesh[LM(i,j+1,k,1)]), 0.50 * (V.Pres[LM(i,j,k,0)] + V.Pres[LM(i,j+1,k,0)]), MESH.Node_Mesh[LM(i,j-1,k,1)], Pressure.Pres[LM(i,j-1,k,0)], MESH.Node_Mesh[LM(i,j,k,1)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j+1,k,1)], Pressure.Pres[LM(i,j+1,k,0)], MESH.Node_Mesh[LM(i,j+2,k,1)], Pressure.Pres[LM(i,j+2,k,0)]);
                
                Pressure.Gradient_Y[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,1)]) * (P_2 - P_1);

                P_1 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k-1,2)] + MESH.Node_Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), MESH.Node_Mesh[LM(i,j,k-2,2)], Pressure.Pres[LM(i,j,k-2,0)], MESH.Node_Mesh[LM(i,j,k-1,2)], Pressure.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], Pressure.Pres[LM(i,j,k+1,0)]);
                P_2 = CS(0.50 * (MESH.Node_Mesh[LM(i,j,k,2)] + MESH.Node_Mesh[LM(i,j,k+1,2)]), 0.50 * (W.Pres[LM(i,j,k,0)] + W.Pres[LM(i,j,k+1,0)]), MESH.Node_Mesh[LM(i,j,k-1,2)], Pressure.Pres[LM(i,j,k-1,0)], MESH.Node_Mesh[LM(i,j,k,2)], Pressure.Pres[LM(i,j,k,0)], MESH.Node_Mesh[LM(i,j,k+1,2)], Pressure.Pres[LM(i,j,k+1,0)], MESH.Node_Mesh[LM(i,j,k+2,2)], Pressure.Pres[LM(i,j,k+2,0)]);

                Pressure.Gradient_Z[LM(i,j,k,0)] = (1.0 / MESH.DeltaP[LM(i,j,k,2)]) * (P_2 - P_1);

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
void CFD_Solver::Get_VelocityGradients(Mesher MESH, Property_Struct &PropertyName){
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

// Function to calculate the dynamic viscosity of the species mix
void CFD_Solver::Get_DynamicViscosity(Species_Solver SPE_S1){
int i, j, k;

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = - 1; k < NZ + 1; k++){
                Stress.mu_Visc[LM(i,j,k,0)] = SPE_S1.JANAF_DynViscosity(T_Pres[LM(i,j,k,0)], i, j, k);
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
                Stress.Tau_xx_X[LMU(i,j,k,0)] =  - (2.0/3.0) * 0.50 * (Stress.mu_Visc[LM(i-1,j,k,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * 0.50 * (Stress.Divergence[LM(i-1,j,k,0)] + Stress.Divergence[LM(i,j,k,0)]) + 2.0 * mu * U.Wall_U_Diff_X[LMU(i,j,k,0)];
            }
        }
    }

    // Tau YY Walls V
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY + 1; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_yy_Y[LMV(i,j,k,0)] =  - (2.0/3.0) * 0.50 * (Stress.mu_Visc[LM(i,j-1,k,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * 0.50 * (Stress.Divergence[LM(i,j-1,k,0)] + Stress.Divergence[LM(i,j,k,0)]) + 2.0 * mu * V.Wall_V_Diff_Y[LMV(i,j,k,0)];
            }
        }
    }

    // Tau ZZ Walls W
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                Stress.Tau_zz_Z[LMW(i,j,k,0)] =  - (2.0/3.0) * 0.50 * (Stress.mu_Visc[LM(i,j,k-1,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * 0.50 * (Stress.Divergence[LM(i,j,k-1,0)] + Stress.Divergence[LM(i,j,k,0)]) + 2.0 * mu * W.Wall_W_Diff_Z[LMW(i,j,k,0)];
            }
        }
    }

    // Tau XY Walls V
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY + 1; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_yx_Y[LMV(i,j,k,0)] =  0.50 * (Stress.mu_Visc[LM(i,j-1,k,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * (U.Wall_V_Diff_Y[LMV(i,j,k,0)] + V.Wall_V_Diff_X[LMV(i,j,k,0)]);
            }
        }
    }

    // Tau ZX Walls W
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                Stress.Tau_zx_Z[LMW(i,j,k,0)] =  0.50 * (Stress.mu_Visc[LM(i,j,k-1,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * (W.Wall_W_Diff_X[LMW(i,j,k,0)] + U.Wall_W_Diff_Z[LMW(i,j,k,0)]);
            }
        }
    }

    // Tau XY Walls U
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_xy_X[LMU(i,j,k,0)] =  0.50 * (Stress.mu_Visc[LM(i-1,j,k,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * (U.Wall_U_Diff_Y[LMU(i,j,k,0)] + V.Wall_U_Diff_X[LMU(i,j,k,0)]);
            }
        }
    }

    // Tau ZY Walls W
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ + 1; k++){
                Stress.Tau_zy_Z[LMW(i,j,k,0)] =  0.50 * (Stress.mu_Visc[LM(i,j,k-1,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * (V.Wall_W_Diff_Z[LMW(i,j,k,0)] + W.Wall_W_Diff_Y[LMW(i,j,k,0)]);
            }
        }
    }

    // Tau YZ Walls V
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY + 1; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_yz_Y[LMV(i,j,k,0)] =  0.50 * (Stress.mu_Visc[LM(i,j-1,k,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * (V.Wall_V_Diff_Z[LMV(i,j,k,0)] + W.Wall_V_Diff_Y[LMV(i,j,k,0)]);
            }
        }
    }

    // Tau XZ Walls U
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Stress.Tau_xz_X[LMU(i,j,k,0)] =  0.50 * (Stress.mu_Visc[LM(i-1,j,k,0)] + Stress.mu_Visc[LM(i,j,k,0)]) * (W.Wall_U_Diff_X[LMU(i,j,k,0)] + U.Wall_U_Diff_Z[LMU(i,j,k,0)]);
            }
        }
    }

}

// Function to calculate the diffusive term of each velocity
void CFD_Solver::Get_DiffusiveTerm_Velocity(Mesher MESH, Property_Struct &PropertyName, double *Wall_U_Term, double *Wall_V_Term, double *Wall_W_Term){
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
void CFD_Solver::Get_StepContribution_Velocity(Property_Struct &PropertyName, double *PressureGradient){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.ContributionPres[LM(i,j,k,0)] = (1.0 / Density.Pres[LM(i,j,k,0)]) * (- PropertyName.Convective[LM(i,j,k,0)] - PressureGradient[LM(i,j,k,0)] + PropertyName.Diffusive[LM(i,j,k,0)]) + PropertyName.Gravity;
            }
        }
    }

}