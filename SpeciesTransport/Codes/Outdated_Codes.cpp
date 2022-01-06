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
