//------------------------------------------------------------------------------------------------//
//                CPP FILE FOR N-S MOMENTUM EQUATION CONVECTION CALCULATIONS                      //
//------------------------------------------------------------------------------------------------//

// Function to calculate the convective term of Velocity U
void CFD_Solver::Get_ConvectionU(Mesher MESH){
int i, j, k;
double uW_pred, uE_pred;
double uW, uE, vS, vN, uS, uN, wH, wT, uH, uT;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                
                uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
                uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

                uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

                vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V.Pres[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V.Pres[LV(i,j,k,0)]);
                vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V.Pres[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V.Pres[LV(i,j+1,k,0)]);

                uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U.Pres[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U.Pres[LU(i,j+2,k,0)], EsquemaLargo);

                wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W.Pres[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W.Pres[LW(i,j,k,0)]);
				wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W.Pres[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W.Pres[LW(i,j,k+1,0)]);

                uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U.Pres[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U.Pres[LU(i,j,k+2,0)], EsquemaLargo);

				U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

            }
        }
    }

    // Top
    j = NY - 1;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (k = 1; k < NZ - 1; k++){
            
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V.Pres[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V.Pres[LV(i,j,k,0)]);
            vN = V.Top[TOP(i,j,k)];

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U.Pres[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], EsquemaLargo);
			uN = U.Top[TOP(i,j,k)];

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W.Pres[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W.Pres[LW(i,j,k,0)]);
			wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W.Pres[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W.Pres[LW(i,j,k+1,0)]);

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U.Pres[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], EsquemaLargo);
			uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U.Pres[LU(i,j,k+2,0)], EsquemaLargo);

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

            
        }
    }

    // Bottom
    j = 0;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (k = 1; k < NZ - 1; k++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

            vS = V.Bottom[BOTTOM(i,j,k)];
            vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V.Pres[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V.Pres[LV(i,j+1,k,0)]);

            uS = U.Bottom[BOTTOM(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U.Pres[LU(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W.Pres[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W.Pres[LW(i,j,k,0)]);
			wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W.Pres[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W.Pres[LW(i,j,k+1,0)]);

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U.Pres[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], EsquemaLargo);
			uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U.Pres[LU(i,j,k+2,0)], EsquemaLargo);

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

        }
    }

    // Here
    k = 0;

    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 1; j < NY - 1; j++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V.Pres[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V.Pres[LV(i,j,k,0)]);
            vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V.Pres[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V.Pres[LV(i,j+1,k,0)]);

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U.Pres[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], EsquemaLargo);
			uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U.Pres[LU(i,j+2,k,0)], EsquemaLargo);

            wH = W.Here[HERE(i,j,k)];
			wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W.Pres[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W.Pres[LW(i,j,k+1,0)]);

            uH = U.Here[HERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U.Pres[LU(i,j,k+2,0)], EsquemaLargo);

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

        }
    }

    // There
    k = NZ - 1;

    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 1; j < NY - 1; j++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V.Pres[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V.Pres[LV(i,j,k,0)]);
            vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V.Pres[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V.Pres[LV(i,j+1,k,0)]);

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U.Pres[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], EsquemaLargo);
			uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U.Pres[LU(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W.Pres[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W.Pres[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U.Pres[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], EsquemaLargo);
			uT = U.There[THERE(i,j,k)];

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

        }
    }

    // Bottom Here Corner
    j = 0;
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
                
        uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
        uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

        uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
		uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

        vS = V.Bottom[BOTTOM(i,j,k)];
        vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V.Pres[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V.Pres[LV(i,j+1,k,0)]);

        uS = U.Bottom[BOTTOM(i,j,k)];
		uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U.Pres[LU(i,j+2,k,0)], EsquemaLargo);

        wH = W.Here[HERE(i,j,k)];
		wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W.Pres[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W.Pres[LW(i,j,k+1,0)]);

        uH = U.Here[HERE(i,j,k)];
		uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U.Pres[LU(i,j,k+2,0)], EsquemaLargo);

		U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }

    // Bottom There Corner
    j = 0;
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

            vS = V.Bottom[BOTTOM(i,j,k)];
            vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V.Pres[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V.Pres[LV(i,j+1,k,0)]);

            uS = U.Bottom[BOTTOM(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U.Pres[LU(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W.Pres[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W.Pres[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U.Pres[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], EsquemaLargo);
			uT = U.There[THERE(i,j,k)];

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }

    // Top Here Corner
    j = NY - 1;
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
                
        uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
        uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

        uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
		uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

        vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V.Pres[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V.Pres[LV(i,j,k,0)]);
        vN = V.Top[TOP(i,j,k)];

        uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U.Pres[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], EsquemaLargo);
		uN = U.Top[TOP(i,j,k)];

        wH = W.Here[HERE(i,j,k)];
		wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W.Pres[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W.Pres[LW(i,j,k+1,0)]);

        uH = U.Here[HERE(i,j,k)];
		uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U.Pres[LU(i,j,k+2,0)], EsquemaLargo);

		U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }

    // Top There Corner
    j = NY - 1;
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U.Pres[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U.Pres[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U.Pres[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V.Pres[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V.Pres[LV(i,j,k,0)]);
            vN = V.Top[TOP(i,j,k)];

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U.Pres[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U.Pres[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U.Pres[LU(i,j+1,k,0)], EsquemaLargo);
			uN = U.Top[TOP(i,j,k)];

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W.Pres[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W.Pres[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U.Pres[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U.Pres[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U.Pres[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U.Pres[LU(i,j,k+1,0)], EsquemaLargo);
			uT = U.There[THERE(i,j,k)];

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }


    if (Rango == 0){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){	
				U.Convective[LU(0,j,k,0)] = 0.0;
			}
		}
    }
    else if (Rango == Procesos - 1){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){	
				U.Convective[LU(NX,j,k,0)] = 0.0;
			}
		}
    }
}