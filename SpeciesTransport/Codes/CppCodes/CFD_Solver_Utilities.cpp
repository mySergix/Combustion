// Function to set initial values for the variables and properties fields
void CFD_Solver::Set_InitialValues(){
int i, j, k;

    // Densities
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Density.Pres[LM(i,j,k,0)] = 0.95;
                Density.Past[LM(i,j,k,0)] = 0.95;
            }
        }
    }

    // Velocities
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){   
            for (k = 0; k < NZ; k++){

                // Velocity U
                U.Pres[LM(i,j,k,0)] = 10.0;
                U.Past[LM(i,j,k,0)] = 10.0;

                // Velocity V
                V.Pres[LM(i,j,k,0)] = 0.0;
                V.Past[LM(i,j,k,0)] = 0.0;

                // Velocity W
                W.Pres[LM(i,j,k,0)] = 0.0;
                W.Past[LM(i,j,k,0)] = 0.0;
            }
        }
    }

    // Pressure
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        for (j = - Halo; j < NY + Halo; j++){
            for (k = - Halo; k < NZ + Halo; k++){
                Pressure.Pres[LM(i,j,k,0)] =  Density.Pres[LM(i,j,k,0)] * 287.0 * 298.15;
            }
        }
    }
    
}

// Function to calulate the time step of the simulation
void CFD_Solver::Get_TimeStep(Mesher MESH){
int i, j, k;
double Courant = 0.2;

DeltaT = 1000.0;

MPI_Status ST;

    // Comparacion
    // Delta += (NewDelta - Delta) * (NewDelta < Delta)

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){

                // CFL Convective

                    // Velocity U
                    DeltaT += (Courant * (MESH.DeltaP[LM(i,j,k,0)] / (abs(U.Pres[LM(i,j,k,0)]) + c + 1e-10)) - DeltaT) * (Courant * (MESH.DeltaP[LM(i,j,k,0)] / (abs(U.Pres[LM(i,j,k,0)]) + c + 1e-10)) < DeltaT);

                    // Velocity V
                    DeltaT += (Courant * (MESH.DeltaP[LM(i,j,k,1)] / (abs(V.Pres[LM(i,j,k,0)]) + c + 1e-10)) - DeltaT) * (Courant * (MESH.DeltaP[LM(i,j,k,1)] / (abs(V.Pres[LM(i,j,k,0)]) + c + 1e-10)) < DeltaT);

                    // Velocity W
                    DeltaT += (Courant * (MESH.DeltaP[LM(i,j,k,2)] / (abs(W.Pres[LM(i,j,k,0)]) + c + 1e-10)) - DeltaT) * (Courant * (MESH.DeltaP[LM(i,j,k,2)] / (abs(W.Pres[LM(i,j,k,0)]) + c + 1e-10)) < DeltaT);
            
            }
        }
    }

    MPI_Allreduce(&DeltaT, &DeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

// Convective Scheme Function
double CFD_Solver::CS(double X, double V, double X1, double Phi1, double X2, double Phi2, double X3, double Phi3, double X4, double Phi4){
double XD, PHID, XC, PHIC, XU, PHIU;
double PhiC_, XC_, XE_, PHIE_, PhiE;

    if (V >= 0){
        XD = X3;
        PHID = Phi3;
        XC = X2;
        PHIC = Phi2;
        XU = X1;
        PHIU = Phi1;
    }
    else{
        XD = X2;
        PHID = Phi2;
        XC = X3;
        PHIC = Phi3;
        XU = X4;
        PHIU = Phi4;
    }

    // Non - Dimensionalization
    PhiC_ = (PHIC - PHIU) / (PHID - PHIU + 1e-10);
    XC_ = (XC - XU) / (XD - XU);
    XE_ = (X - XU) / (XD - XU);

    // CDS Scheme
    PHIE_ = ((XE_ - XC_) / (1 - XC_)) + ((XE_ - 1) / (XC_ - 1)) * PhiC_;

    // Dimensionalization
    PhiE = PHIU + (PHID - PHIU) * PHIE_;

    return PhiE;

}

// Calculation of Property Value at the Walls of the Volumes
void CFD_Solver::Get_WallsValue_Scalar(Parallel P1, double *Mesh, Property &PropertyName){
int i, j, k;

    // Communication of the local metrix of the properties
    P1.CommunicateLocalMatrix(PropertyName.Pres, PropertyName.Pres);

    // MU Walls
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Wall_U[LMU(i,j,k,0)] = CS(0.50 * (Mesh[LM(i - 1,j,k,0)] + Mesh[LM(i,j,k,0)]), 0.50 * (U.Pres[LM(i - 1,j,k,0)] + U.Pres[LM(i,j,k,0)]), Mesh[LM(i-2,j,k,0)], PropertyName.Pres[LM(i-2,j,k,0)], Mesh[LM(i-1,j,k,0)], PropertyName.Pres[LM(i-1,j,k,0)], Mesh[LM(i,j,k,0)], PropertyName.Pres[LM(i,j,k,0)], Mesh[LM(i+1,j,k,0)], PropertyName.Pres[LM(i+1,j,k,0)]); 
            }
        }
    }

    // MV Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Wall_V[LMV(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j-1,k,1)] + Mesh[LM(i,j,k,1)]), 0.50 * (V.Pres[LM(i,j-1,k,0)] + V.Pres[LM(i,j,k,0)]), Mesh[LM(i,j-2,k,1)], PropertyName.Pres[LM(i,j-2,k,0)], Mesh[LM(i,j-1,k,1)], PropertyName.Pres[LM(i,j-1,k,0)], Mesh[LM(i,j,k,1)], PropertyName.Pres[LM(i,j,k,0)], Mesh[LM(i,j+1,k,1)], PropertyName.Pres[LM(i,j+1,k,0)]);       
            }
        }
    }

    // MW Walls
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 1; k < NZ; k++){
                PropertyName.Wall_W[LMW(i,j,k,0)] = CS(0.50 * (Mesh[LM(i,j,k-1,2)] + Mesh[LM(i,j,k,2)]), 0.50 * (W.Pres[LM(i,j,k-1,0)] + W.Pres[LM(i,j,k,0)]), Mesh[LM(i,j,k-2,2)], PropertyName.Pres[LM(i,j,k-2,0)], Mesh[LM(i,j,k-1,2)], PropertyName.Pres[LM(i,j,k-1,0)], Mesh[LM(i,j,k,2)], PropertyName.Pres[LM(i,j,k,0)], Mesh[LM(i,j,k+1,2)], PropertyName.Pres[LM(i,j,k+1,0)]);
            }
        }
    }

    ApplyBoundaries(PropertyName);   

}

// Function to integrate the density calculation
void CFD_Solver::Get_TemporalIntegration(Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Fut[LM(i,j,k,0)] = (2.0 * Beta * PropertyName.Pres[LM(i,j,k,0)] - (Beta - 0.50) * PropertyName.Past[LM(i,j,k,0)]) / (Beta + 0.50) 
                                              + DeltaT * ((1.0 + Beta) * PropertyName.ContributionPres[LM(i,j,k,0)] - Beta * PropertyName.ContributionPast[LM(i,j,k,0)]); 
            }
        }
    }

}

// Function to update the values of the properties fields
void CFD_Solver::UpdateField(Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.Past[LM(i,j,k,0)] = PropertyName.Pres[LM(i,j,k,0)];
                PropertyName.Pres[LM(i,j,k,0)] = PropertyName.Fut[LM(i,j,k,0)];

                PropertyName.ContributionPast[LM(i,j,k,0)] = PropertyName.ContributionPres[LM(i,j,k,0)];
            }
        }
    }

}

// Function to calculate the difference between time steps
void CFD_Solver::Get_MaximumDifference(Property &PropertyName, double &MaxDiff){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                MaxDiff += (abs((PropertyName.Fut[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j,k,0)]) / (PropertyName.Pres[LM(i,j,k,0)] + 1e-10)) - MaxDiff) * (abs((PropertyName.Fut[LM(i,j,k,0)] - PropertyName.Pres[LM(i,j,k,0)]) / (PropertyName.Pres[LM(i,j,k,0)] + 1e-10)) > MaxDiff);
            }
        }
    }

    MPI_Allreduce(&MaxDiff, &MaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

}

// Function to Check for convergence criteria
void CFD_Solver::Get_ConvergenceCriteria(){
int i, j, k;
MaxDiff = 0.0;

    // Check Convergence for Density
    Get_MaximumDifference(Density, MaxDiff);

}
