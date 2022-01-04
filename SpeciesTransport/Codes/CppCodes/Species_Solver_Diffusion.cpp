//------------------------------------------------------------------------------------------------//
//                          CPP FILE FOR SPECIES DIFFUSION CALCULATIONS                           //
//------------------------------------------------------------------------------------------------//

// Function to calculate the binary diffusion coefficient based on Champan-Enskog model
double Species_Solver::Get_BinaryDiff_ChampanEnskog(int SP, double Temperature, double Pressure, int i, int j, int k){
int i;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn;

    for (i = 0; i < N_Species; i++){
        if (i =! SP){
            W_ab = power(1.0 / Species[SP].W_molar + 1.0 / Species[i].W_molar, - 1.0);
            Sigma_ab = (Species[SP].sigma + Species[i].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[SP].Epsilon * Species[i].Epsilon);
            E_ab = Epsilon_ab / kB;
            Tn = Temperature / E_ab;
            OmegaD = 1.06036 / power(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_ab = 10.1325 * (0.001858 * power(Tn, 1.5) * power(W_ab, -0.5)) / (Pressure * power(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[i].X[LM(i,j,k,0)] / D_ab;
        }    
    }
       
    return Sum;

}

// Function to calculate the binary difussion coefficient based on Wilke-Lee model 
double Species_Solver::Get_BinaryDiff_WilkeLee(int SP, double Temperature, double Pressure, int i, int j, int k){
int i;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn;

    for (i = 0; i < N_Species; i++){
        if (i =! SP){
            W_ab = power(1.0 / Species[SP].W_molar + 1.0 / Species[i].W_molar, - 1.0);
            Sigma_ab = (Species[SP].sigma + Species[i].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[SP].Epsilon * Species[i].Epsilon);
            E_ab = Epsilon_ab / kB;
            Tn = Temperature / E_ab;
            OmegaD = 1.06036 / power(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_ab = 10.1325 * ((0.0027 - 0.0005 * power(W_ab, - 0.5)) * power(Tn, 1.5) * power(W_ab, -0.5)) / (Pressure * power(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[i].X[LM(i,j,k,0)] / D_ab;
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
                Species[SP].D_am[LM(i,j,k,0)] = (1.0 - Species[SP].X[LM(i,j,k,0)]) / Get_BinaryDif_ChampanEnskog(SP, CFD_S1.T.Pres[LP(i,j,k,0)], CFD_S1.Pressure.Pres[LP(i,j,k,0)], i, j, k);
            }
        }
    }

}

// Function to calculate the walls value of diffusion coefficients
void Species_Solver::Get_WallsDiffusionCoefficients(int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].D_am_Wall_U[LMU(i,j,k,0)] = 0.50 * (Species[SP].D_am[LM(i-1,j,k,0)] + Species[SP].D_am[LM(i,j,k,0)]);
                // Faltan direcciones Y y Z
            }
        }
    }

}

// Function to calculate the diffusion term of a species
void Species_Solver::Get_SpeciesDiffusion(Mesher MESH, CFD_Solver CFD_S1, int SP){
int i, j, k;

    // Cambiarlo de tal forma que se calcule la velocidad de difusion
    // Aqui esta calculada como D * grad(y)
    // Ponerla para calcularla en funcion de X (fraccion molar) (Poinsot)
    // Importante, la aproximacion solo vale para flujos sin fuerzas volumetricas y SIN GRADIENTES DE PRESION
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].DiffusionTerm[LM(i,j,k,0)] = (1.0 / MESH.VolMP) * (
                                                       - MESH.Surf[LM(i,j,k,0)] * CFD_S1.Density.Wall_U[LMU(i,j,k,0)] * Species[SP].D_am_Wall_U[LMU(i,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i-1,j,k,0)] + MESH.DeltaP[LM(i,j,k,0)])) * (Species[SP].Y[LM(i,j,k,0)] - Species[SP].Y[LM(i-1,j,k,0)])
                                                       + MESH.Surf[LM(i,j,k,0)] * CFD_S1.Density.Wall_U[LMU(i+1,j,k,0)] * Species[SP].D_am_Wall_U[LMU(i+1,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k,0)] + MESH.DeltaP[LM(i+1,j,k,0)])) * (Species[SP].Y[LM(i+1,j,k,0)] - Species[SP].Y[LM(i,j,k,0)])
                                                       - MESH.Surf[LM(i,j,k,1)] * CFD_S1.Density.Wall_V[LMV(i,j,k,0)] * Species[SP].D_am_Wall_V[LMV(i,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j-1,k,1)] + MESH.DeltaP[LM(i,j,k,1)])) * (Species[SP].Y[LM(i,j,k,0)] - Species[SP].Y[LM(i,j-1,k,0)])
                                                       + MESH.Surf[LM(i,j,k,1)] * CFD_S1.Density.Wall_V[LMV(i,j+1,k,0)] * Species[SP].D_am_Wall_V[LMV(i,j+1,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k,1)] + MESH.DeltaP[LM(i,j+1,k,1)])) * (Species[SP].Y[LM(i,j+1,k,0)] - Species[SP].Y[LM(i,j,k,0)])
                                                       - MESH.Surf[LM(i,j,k,2)] * CFD_S1.Density.Wall_W[LMW(i,j,k,0)] * Species[SP].D_am_Wall_W[LMW(i,j,k,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k-1,2)] + MESH.DeltaP[LM(i,j,k,2)])) * (Species[SP].Y[LM(i,j,k,0)] - Species[SP].Y[LM(i,j,k-1,0)])
                                                       + MESH.Surf[LM(i,j,k,2)] * CFD_S1.Density.Wall_W[LMW(i,j,k+1,0)] * Species[SP].D_am_Wall_W[LMW(i,j,k+1,0)] * (2.0 / (MESH.DeltaP[LM(i,j,k,2)] + MESH.DeltaP[LM(i,j,k+1,2)])) * (Species[SP].Y[LM(i,j,k+1,0)] - Species[SP].Y[LM(i,j,k,0)])
                                                       );
            }
        }
    }

}