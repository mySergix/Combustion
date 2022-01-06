//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR N-S MASS EQUATION CALCULATIONS                        //
//------------------------------------------------------------------------------------------------//

// Function to compute the Convective Term of Density
void CFD_Solver::Get_ConvectiveTerm_Density(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Density.ContributionPres[LM(i,j,k,0)] = - (1.0/MESH.Vol[LM(i,j,k,0)]) * (
                                                     + MESH.Surf[LM(i,j,k,0)] * Density.Wall_U[LMU(i+1,j,k,0)] * U.Wall_U[LMU(i+1,j,k,0)]
                                                     - MESH.Surf[LM(i,j,k,0)] * Density.Wall_U[LMU(i,j,k,0)] * U.Wall_U[LMU(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,1)] * Density.Wall_V[LMV(i,j+1,k,0)] * V.Wall_V[LMV(i,j+1,k,0)]
                                                     - MESH.Surf[LM(i,j,k,1)] * Density.Wall_V[LMV(i,j,k,0)] * V.Wall_V[LMV(i,j,k,0)]
                                                     + MESH.Surf[LM(i,j,k,2)] * Density.Wall_W[LMW(i,j,k+1,0)] * W.Wall_W[LMW(i,j,k+1,0)]
                                                     - MESH.Surf[LM(i,j,k,2)] * Density.Wall_W[LMW(i,j,k,0)] * W.Wall_W[LMW(i,j,k,0)]
                                                     );
            }
        }
    }

}
