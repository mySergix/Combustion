//------------------------------------------------------------------------------------------------//
//                         CPP FILE FOR N-S ENERGY EQUATION CALCULATIONS                          //
//------------------------------------------------------------------------------------------------//

// Function to calculate the initial absolute enthalpy of the mix at the beginning
void CFD_Solver::Get_Initial_AbsEnthalpy(Species_Solver SPE_S1){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Hs.Pres[LM(i,j,k,0)] = SPE_S1.JANAF_AbsEnthalpy_Specie_Mix(T.Pres[LM(i,j,k,0)], i, j, k);
                Hs.Past[LM(i,j,k,0)] = Hs.Pres[LM(i,j,k,0)];
            }
        }
    }

}

// Function to calculate the step contribution to each velocity
void CFD_Solver::Get_StepContribution_Enthalpy(){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Hs.ContributionPres[LM(i,j,k,0)] = (1.0 / Density.Pres[LM(i,j,k,0)]) * (- Hs.Convective[LM(i,j,k,0)]);
            }
        }
    }

}

// Function to calculate the thermal conductivity of the domain
void CFD_Solver::Get_ThermalConductivity(Species_Solver SPE_S1){
int i, j, k;

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = - 1; k < NZ + 1; k++){ // Vigilar -1 y +1 en la condiciones de contorno
                Lambda[LM(i,j,k,0)] = SPE_S1.JANAF_ThermalCond(T.Pres[LM(i,j,k,0)], i, j, k);
            }
        }
    }

}

// Function to calculate the Fourier term of the energy equation
void CFD_Solver::Get_FourierDiffusion(Mesher MESH){
int i, j, k;

    // Falta aplicar las condiciones de contorno en la Z y la Y
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Hs.FourierDiffusion[LM(i,j,k,0)] = (1.0 / MESH.Vol[LM(i,j,k,0)]) * (
                                                 - MESH.Surf[LM(i,j,k,0)] * 0.50 * (Lambda[LM(i-1,j,k,0)] + Lambda[LM(i,j,k,0)]) * (2.0 / (MESH.DeltaP[LM(i-1,j,k,0)] + MESH.DeltaP[LM(i,j,k,0)])) * (T.Pres[LM(i,j,k,0)] - T.Pres[LM(i-1,j,k,0)])
                                                 + MESH.Surf[LM(i,j,k,0)] * 0.50 * (Lambda[LM(i,j,k,0)] + Lambda[LM(i+1,j,k,0)]) * (2.0 / (MESH.DeltaP[LM(i,j,k,0)] + MESH.DeltaP[LM(i+1,j,k,0)])) * (T.Pres[LM(i+1,j,k,0)] - T.Pres[LM(i,j,k,0)])
                                                 - MESH.Surf[LM(i,j,k,1)] * 0.50 * (Lambda[LM(i,j-1,k,0)] + Lambda[LM(i,j,k,0)]) * (2.0 / (MESH.DeltaP[LM(i,j-1,k,1)] + MESH.DeltaP[LM(i,j,k,1)])) * (T.Pres[LM(i,j,k,0)] - T.Pres[LM(i,j-1,k,0)])
                                                 + MESH.Surf[LM(i,j,k,1)] * 0.50 * (Lambda[LM(i,j,k,0)] + Lambda[LM(i,j+1,k,0)]) * (2.0 / (MESH.DeltaP[LM(i,j,k,1)] + MESH.DeltaP[LM(i,j+1,k,1)])) * (T.Pres[LM(i,j+1,k,0)] - T.Pres[LM(i,j,k,0)])
                                                 - MESH.Surf[LM(i,j,k,2)] * 0.50 * (Lambda[LM(i-1,j,k,0)] + Lambda[LM(i,j,k,0)]) * (2.0 / (MESH.DeltaP[LM(i,j,k-1,2)] + MESH.DeltaP[LM(i,j,k,2)])) * (T.Pres[LM(i,j,k,0)] - T.Pres[LM(i,j,k-1,0)])
                                                 + MESH.Surf[LM(i,j,k,2)] * 0.50 * (Lambda[LM(i,j,k,0)] + Lambda[LM(i,j,k+1,0)]) * (2.0 / (MESH.DeltaP[LM(i,j,k,2)] + MESH.DeltaP[LM(i,j,k+1,2)])) * (T.Pres[LM(i,j,k+1,0)] - T.Pres[LM(i,j,k,0)])
                                                 );
            }
        }
    }

}