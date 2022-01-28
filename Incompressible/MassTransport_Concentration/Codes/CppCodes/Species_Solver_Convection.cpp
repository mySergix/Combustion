//------------------------------------------------------------------------------------------------//
//                   CPP FILE FOR SPECIES EQUATION CONVECTION CALCULATIONS                        //
//------------------------------------------------------------------------------------------------//

// Function to calculate the convective scheme
inline double Species_Solver::ConvectiveSchemeSpecies(double CoordObjetivo, double Velocity, double Coord1, double Phi1, double Coord2, double Phi2, double Coord3, double Phi3, double Coord4, double Phi4, string Esquema){

double PhiObjetivo;

double CoordD;
double PhiD;

double CoordC;
double PhiC;

double CoordU;
double PhiU;

	if (Velocity <= 0.0 || (Phi1 == 0.0 && Coord1 == 0.0)){

		CoordD = Coord2;
		PhiD = Phi2;
		CoordC = Coord3;
		PhiC = Phi3;
		CoordU = Coord4;
		PhiU = Phi4;

	}
	else if(Velocity > 0.0 || (Phi4 == 0.0 && Coord4 == 0.0)){

		CoordD = Coord3;
		PhiD = Phi3;
		CoordC = Coord2;
		PhiC = Phi2;
		CoordU = Coord1;
		PhiU = Phi1;

	}

	//Adimensionalizacion
	double PhiAdimC = (PhiC - PhiU)/(PhiD - PhiU);
	double AdimCoordC = (CoordC - CoordU)/(CoordD - CoordU);
	double AdimCoordE = (CoordObjetivo - CoordU)/(CoordD - CoordU);

	//Evaluacion
	double PhiF;

	if (PhiD == PhiU){
		PhiObjetivo = PhiD;
	}
	else{

        // CDS
        PhiF = ((AdimCoordE - AdimCoordC) / (1 - AdimCoordC)) + ((AdimCoordE - 1) / (AdimCoordC - 1)) * PhiAdimC;

        // QUICK
		//PhiF = AdimCoordE + (((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0))))*(PhiAdimC - AdimCoordC);

		//Dimensionalizacion
		PhiObjetivo = PhiU + (PhiD - PhiU)*PhiF;
	}

	return PhiObjetivo;

}

// Function to calculate the convection of each species
void Species_Solver::Get_ConvectionSpecies(Mesher MESH, CFD_Solver CFD_S1, int SP){
int i, j, k;
double Cw, Ce, Cs, Cn, Ch, Ct;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

            }
        }
    }

    // Bottom
    j = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 1; k < NZ - 1; k++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = Species[SP].Bottom[BOTTOM(i,j,k)];
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }
    }

    // Top
    j = NY - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 1; k < NZ - 1; k++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = Species[SP].Top[TOP(i,j,k)];

            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }
    }

    // Here
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY - 1; j++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = Species[SP].Here[HERE(i,j,k)];
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }
    }

    // There
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY - 1; j++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = Species[SP].There[THERE(i,j,k)];

            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }
    }

    if (Rango == 0){
        i = 0;

        // Center
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                Cw = Species[SP].Left[LEFT(i,j,k)];
                Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

            }
        }

        // Bottom
        j = 0;
        for (k = 1; k < NZ - 1; k++){
            Cw = Species[SP].Left[LEFT(i,j,k)];
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = Species[SP].Bottom[BOTTOM(i,j,k)];
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){
            Cw = Species[SP].Left[LEFT(i,j,k)];
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = Species[SP].Top[TOP(i,j,k)];

            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

        // Here
        k = 0;
        for (j = 1; j < NY - 1; j++){
            Cw = Species[SP].Left[LEFT(i,j,k)];
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = Species[SP].Here[HERE(i,j,k)];
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

        // There
        k = NZ - 1;
        for (j = 1; j < NY - 1; j++){
            Cw = Species[SP].Left[LEFT(i,j,k)];
            Ce = ConvectiveSchemeSpecies(MESH.MU[LU(i+1,j,k,0)], CFD_S1.U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[SP].C_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = Species[SP].There[THERE(i,j,k)];

            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

    }
    else if (Rango == Procesos - 1){
        i = NX - 1;

        // Center
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ce = Species[SP].Right[RIGHT(i,j,k)];
                
                Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

            }
        }

        // Bottom
        j = 0;
        for (k = 1; k < NZ - 1; k++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = Species[SP].Right[RIGHT(i,j,k)];

            Cs = Species[SP].Bottom[BOTTOM(i,j,k)];
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = Species[SP].Right[RIGHT(i,j,k)];

            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = Species[SP].Top[TOP(i,j,k)];

            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

        // Here
        k = 0;
        for (j = 1; j < NY - 1; j++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = Species[SP].Right[RIGHT(i,j,k)];

            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = Species[SP].Here[HERE(i,j,k)];
            Ct = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k+1,2)], CFD_S1.W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[SP].C_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

        // There
        k = NZ - 1;
        for (j = 1; j < NY - 1; j++){
            Cw = ConvectiveSchemeSpecies(MESH.MU[LU(i,j,k,0)], CFD_S1.U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[SP].C_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[SP].C_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[SP].C_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ce = Species[SP].Right[RIGHT(i,j,k)];

            Cs = ConvectiveSchemeSpecies(MESH.MV[LV(i,j,k,1)], CFD_S1.V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[SP].C_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Cn = ConvectiveSchemeSpecies(MESH.MV[LV(i,j+1,k,1)], CFD_S1.V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[SP].C_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[SP].C_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[SP].C_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Ch = ConvectiveSchemeSpecies(MESH.MW[LW(i,j,k,2)], CFD_S1.W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[SP].C_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[SP].C_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[SP].C_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[SP].C_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Ct = Species[SP].There[THERE(i,j,k)];

            Species[SP].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (CFD_S1.U.Pres[LU(i+1,j,k,0)] * Ce - Cw * CFD_S1.U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (CFD_S1.V.Pres[LV(i,j+1,k,0)] * Cn - Cs * CFD_S1.V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (CFD_S1.W.Pres[LW(i,j,k+1,0)] * Ct - Ch * CFD_S1.W.Pres[LW(i,j,k,0)])
											 );

        }

    }

}

