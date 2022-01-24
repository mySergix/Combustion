//------------------------------------------------------------------------------------------------//
//                          CPP FILE FOR SPECIES DIFFUSION CALCULATIONS                           //
//------------------------------------------------------------------------------------------------//

// Function to calculate the binary diffusion coefficient based on Champan-Enskog model
double Species_Solver::Get_BinaryDiff_ChampanEnskog(int SP, double Temperature, double Pressure, int i, int j, int k){
int n;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn, OmegaD;
double D_ab;

    for (n = 0; n < N_Species; n++){
        if (n =! SP){
            W_ab = pow(1.0 / Species[SP].Wmolar + 1.0 / Species[n].Wmolar, - 1.0);
            Sigma_ab = (Species[SP].sigma + Species[n].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[SP].Epsilon * Species[n].Epsilon);
            Tn = Temperature / Epsilon_ab;
            OmegaD = 1.06036 / pow(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_ab = 10.1325 * (0.001858 * pow(Tn, 1.5) * pow(W_ab, -0.5)) / (Pressure * pow(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[n].X[LM(i,j,k,0)] / D_ab;
        }    
    }
       
    return Sum;

}

// Function to calculate the binary difussion coefficient based on Wilke-Lee model 
double Species_Solver::Get_BinaryDiff_WilkeLee(int SP, double Temperature, double Pressure, int i, int j, int k){
int n;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn, OmegaD;
double D_ab;

    for (n = 0; n < N_Species; n++){
        if (n =! SP){
            W_ab = pow(1.0 / Species[SP].Wmolar + 1.0 / Species[n].Wmolar, - 1.0);
            Sigma_ab = (Species[SP].sigma + Species[n].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[SP].Epsilon * Species[n].Epsilon);
            Tn = Temperature / Epsilon_ab;
            OmegaD = 1.06036 / pow(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_ab = 10.1325 * ((0.0027 - 0.0005 * pow(W_ab, - 0.5)) * pow(Tn, 1.5) * pow(W_ab, -0.5)) / (Pressure * pow(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[n].X[LM(i,j,k,0)] / D_ab;
        }    
    }
       
    return Sum;

}

// Function to calculate the diffusion coefficient of a species
void Species_Solver::Get_DiffusionCoefficient_FickModel(int SP, CFD_Solver CFD_S1){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].D_am[LM(i,j,k,0)] = (1.0 - Species[SP].Y_Pres[LM(i,j,k,0)]) / Get_BinaryDiff_ChampanEnskog(SP, CFD_S1.T_Pres[LM(i,j,k,0)], CFD_S1.Pressure.Pres[LM(i,j,k,0)], i, j, k);
            }
        }
    }

}

// Function to calculate the diffusion velocities of the species
void Species_Solver::Get_WallsDiffusionVelocities(int SP, Mesher MESH){
int i, j, k;

    // Hirschfelder and Curtiss aproximation

    // U Diffusion Velocity
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].U_Diff[LMU(i,j,k,0)] = - ((0.50 * (Species[SP].D_am[LM(i-1,j,k,0)] + Species[SP].D_am[LM(i,j,k,0)])) / (0.50 * (Species[SP].X[LM(i-1,j,k,0)] + Species[SP].X[LM(i,j,k,0)]))) * ((Species[SP].X[LM(i,j,k,0)] - Species[SP].X[LM(i-1,j,k,0)]) / (2.0 / (MESH.DeltaP[LM(i-1,j,k,0)] + MESH.DeltaP[LM(i,j,k,0)])));
            }
        }
    }

    // V Diffusion Velocity
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY - 1; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].V_Diff[LMV(i,j,k,0)] = - ((0.50 * (Species[SP].D_am[LM(i,j-1,k,0)] + Species[SP].D_am[LM(i,j,k,0)])) / (0.50 * (Species[SP].X[LM(i,j-1,k,0)] + Species[SP].X[LM(i,j,k,0)]))) * ((Species[SP].X[LM(i,j,k,0)] - Species[SP].X[LM(i,j-1,k,0)]) / (2.0 / (MESH.DeltaP[LM(i,j-1,k,1)] + MESH.DeltaP[LM(i,j,k,1)])));
            }
        }
    }

    // W Diffusion Velocity
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].W_Diff[LMW(i,j,k,0)] = - ((0.50 * (Species[SP].D_am[LM(i,j,k-1,0)] + Species[SP].D_am[LM(i,j,k,0)])) / (0.50 * (Species[SP].X[LM(i,j,k,0)] + Species[SP].X[LM(i,j,k-1,0)]))) * ((Species[SP].X[LM(i,j,k,0)] - Species[SP].X[LM(i,j,k-1,0)]) / (2.0 / (MESH.DeltaP[LM(i,j,k-1,2)] + MESH.DeltaP[LM(i,j,k,2)])));
            }
        }
    }

}