//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES DIFFUSION CALCULATIONS                              //
//------------------------------------------------------------------------------------------------//

// Function to calculate species binary diffusion coefficients
void Species_Solver::Get_DiffusionCoefficients(int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].D_ab[LP(i,j,k,0)] = D_AB;
            }
        }
    }

}

// Function to calculate the species diffusion
void Species_Solver::Get_DiffusionSpecies(Mesher MESH, int SP){
int i, j, k;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }
        }
    }

    // Bottom
    j = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 1; k < NZ - 1; k++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }
    }

    // Top
    j = NY - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 1; k < NZ - 1; k++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].Top[TOP(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }
    }

    // Here
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY - 1; j++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }
    }

    // There
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY - 1; j++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].There[THERE(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }
    }

    if (Rango == 0){
        i = 0;

        // Center
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }
        }
    
        // Bottom
        j = 0;
        for (k = 1; k < NZ - 1; k++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].Top[TOP(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // Here
        k = 0;
        for (j = 1; j < NY - 1; j++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // There
        k = NZ - 1;
        for (j = 1; j < NY - 1; j++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i+1,j,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].There[THERE(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

    }
    else if (Rango == Procesos - 1){
        i = NX - 1;

        // Center
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].Right[RIGHT(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }
        }

        // Bottom
        j = 0;
        for (k = 1; k < NZ - 1; k++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].Right[RIGHT(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].Right[RIGHT(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].Top[TOP(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // Here
        k = 0;
        for (j = 1; j < NY - 1; j++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].Right[RIGHT(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k+1,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // There
        k = NZ - 1;
        for (j = 1; j < NY - 1; j++){
            Species[SP].Diffusive[LP(i,j,k,0)] = Species[SP].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[SP].Right[RIGHT(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j+1,k,0)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[SP].There[THERE(i,j,k)] - Species[SP].C_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[SP].C_Pres[LP(i,j,k,0)] - Species[SP].C_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

    }

}

// Function to set the last species concentration to be (1 - sum(C))
void Species_Solver::Get_MassConservationSpecies(int SP){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[SP].C_Fut[LP(i,j,k,0)] = 1.0 - Species[0].C_Fut[LP(i,j,k,0)];
            }
        }
    }

}