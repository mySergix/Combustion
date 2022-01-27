//------------------------------------------------------------------------------------------------//
//                   CPP FILE FOR CFD SOLVER CLASS POISSON COEFFICIENTS                           //
//------------------------------------------------------------------------------------------------//

// Function to calculate the Poisson Coefficients
void CFD_Solver::Get_PoissonCoefficients(Mesher MESH){
int i, j, k;

    // West and East Coefficients (aw, ae)
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){	
                A.aw[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,0)]/MESH.DeltasMU[LU(i,j,k,0)];
                A.ae[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,0)]/MESH.DeltasMU[LU(i+1,j,k,0)];
            }
        }
    }

    if (Rango == 0){
        for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
                A.aw[LA(0,j,k,0)] = 0.0;
            }
        }
    }
    else if (Rango == Procesos - 1){
        for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
                A.ae[LA(NX-1,j,k,0)] = 0.0;
            }
        }
    }

    // South and North Coefficients (as, an)
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 0; k < NZ; k++){
            
            // Top
            A.as[LA(i,NY-1,k,0)] = MESH.SupMP[LP(i,NY-1,k,1)]/MESH.DeltasMV[LV(i,NY-1,k,1)];
            A.an[LA(i,NY-1,k,0)] = 0.0;

            // Bottom
            A.as[LA(i,0,k,0)] = 0.0;
            A.an[LA(i,0,k,0)] = MESH.SupMP[LP(i,0,k,1)]/MESH.DeltasMV[LV(i,1,k,1)];

			for(j = 1; j < NY - 1; j++){
                A.as[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,1)]/MESH.DeltasMV[LV(i,j,k,1)];
                A.an[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,1)]/MESH.DeltasMV[LV(i,j+1,k,1)];
            }
        }
    }

    // Here and There Coefficients (ah, at)
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){

            // Here
            A.ah[LA(i,j,0,0)] = 0.0;
			A.at[LA(i,j,0,0)] = MESH.SupMP[LP(i,j,0,2)]/MESH.DeltasMW[LW(i,j,1,2)];

            // There
            A.ah[LA(i,j,NZ-1,0)] = MESH.SupMP[LP(i,j,NZ-1,2)]/MESH.DeltasMW[LW(i,j,NZ-1,2)];
			A.at[LA(i,j,NZ-1,0)] = 0.0;

            for(k = 1; k < NZ - 1; k++){
                A.ah[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,2)]/MESH.DeltasMW[LW(i,j,k,2)];
				A.at[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,2)]/MESH.DeltasMW[LW(i,j,k+1,2)];
            }
        }
    }

	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
                A.ap[LA(i,j,k,0)] = A.aw[LA(i,j,k,0)] + A.ae[LA(i,j,k,0)] + A.as[LA(i,j,k,0)] + A.an[LA(i,j,k,0)] + A.ah[LA(i,j,k,0)] + A.at[LA(i,j,k,0)];
            }
        }
    }

}

// Function to write a .txt with the mesh data
void CFD_Solver::PrintTxt(){
int i, j, k;
string Carpeta = "GnuPlotResults/";
ofstream file;
//string FileName;
stringstream InitialNameMP;
string FinalNameMP;
	char FileName[300]; 
	sprintf(FileName, "ap_Processor_%d.txt", Rango);
	//FileName = "DeltasMP_X.txt";

	InitialNameMP<<"../"<<Carpeta<<FileName;

	FinalNameMP = InitialNameMP.str();
    file.open(FinalNameMP.c_str());
	file<<"Processor: "<<Rango<<endl;
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				file<<"I: "<<i<<"\t"<<"J: "<<j<<"\t"<<"K: "<<k<<"\t"<<"ap: "<<A.ap[LA(i,j,k,0)]<<"\t"<<endl;
			}
			file<<endl;
		}   	
    }

	file.close();

}


// Function to solve the system with Gauss-Seidel
void CFD_Solver::Get_GaussSeidel(Parallel P1){
int i, j, k;
MaxDiffGS = 2.0*ConvergenciaGS;

    while (MaxDiffGS >= ConvergenciaGS){

        // Pressure Halo Communication
        P1.CommunicateDataLP(P.Pres, P.Pres);

        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 1; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    P.Pres[LP(i,j,k,0)] = (A.aw[LA(i,j,k,0)]*P.Pres[LP(i-1,j,k,0)] + A.ae[LA(i,j,k,0)]*P.Pres[LP(i+1,j,k,0)] + A.as[LA(i,j,k,0)]*P.Pres[LP(i,j-1,k,0)] + A.an[LA(i,j,k,0)]*P.Pres[LP(i,j+1,k,0)] + A.ah[LA(i,j,k,0)]*P.Pres[LP(i,j,k-1,0)] + A.at[LA(i,j,k,0)]*P.Pres[LP(i,j,k+1,0)] + A.bp[LA(i,j,k,0)])/A.ap[LA(i,j,k,0)];
                }
            }
        }

        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (k = 0; k < NZ; k++){
                P.Pres[LP(i,0,k,0)] = 0.0;
            }
        }

        MaxDiffGS = 0.0;

		for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = 0; j < NY; j++){
			    for(k = 0; k < NZ; k++){	
                    MaxDiffGS += (abs(P.Pres[LP(i,j,k,0)] - P.Sup[LA(i,j,k,0)]) - MaxDiffGS) * (abs(P.Pres[LP(i,j,k,0)] - P.Sup[LA(i,j,k,0)]) >= MaxDiffGS);
					P.Sup[LA(i,j,k,0)] = P.Pres[LP(i,j,k,0)];
				}
			}
		}

        MPI_Allreduce(&MaxDiffGS, &MaxDiffGS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
    }

}