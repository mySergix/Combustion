#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>
//#include <chrono>

using namespace std;

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/ParPro.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/Solver.h"

#define PI 3.141592653589793



//Condiciones de Contorno

//Parte Izquierda
#define PLEFT(i,j,k) NY*(k) + (j)
#define ULEFT(i,j,k) NY*(k) + (j) 
#define VLEFT(i,j,k) (NY+1)*(k) + (j)
#define WLEFT(i,j,k) NY*(k) + (j)

//Parte Derecha
#define PRIGHT(i,j,k) NY*(k) + (j)
#define URIGHT(i,j,k) NY*(k) + (j)
#define VRIGHT(i,j,k) (NY+1)*(k) + (j)
#define WRIGHT(i,j,k) NY*(k) + (j)

//Parte Inferior
#define PBOT(i,j,k) NZ*(i - Ix) + k
#define UBOT(i,j,k) NZ*(i - Ix) + k
#define VBOT(i,j,k) NZ*(i - Ix + 1) + k
#define WBOT(i,j,k) (NZ+1)*(i - Ix + 1) + k

//Parte Superior
#define PTOP(i,j,k) NZ*(i - Ix) + k
#define UTOP(i,j,k) NZ*(i - Ix) + k
#define VTOP(i,j,k) NZ*(i - Ix + 1) + k
#define WTOP(i,j,k) (NZ+1)*(i - Ix + 1) + k

//Parte Here
#define PHERE(i,j,k) NY*(i - Ix) + j
#define UHERE(i,j,k) NY*(i - Ix) + j
#define VHERE(i,j,k) (NY+1)*(i - Ix + 1) + j
#define WHERE(i,j,k) NY*(i - Ix + 1) + j
 
//Parte There 
#define PTHERE(i,j,k) NY*(i - Ix) + j 
#define UTHERE(i,j,k) NY*(i - Ix) + j
#define VTHERE(i,j,k) (NY+1)*(i - Ix + 1) + j
#define WTHERE(i,j,k) NY*(i - Ix + 1) + j

#define Interpolacion(CoordObjetivo, Coord1, Valor1, Coord2, Valor2) Valor1 + ((Valor2 - Valor1)/(Coord2 - Coord1))*(CoordObjetivo - Coord1) //Interpolación Lineal

//Constructor del mallador
Solver::Solver(Memory M1, ReadData R1, ParPro MPI1, Mesher MESH, PostProcessing POST1, string InputDirectorio){
		
	
	//Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	NX = MESH.NX;
	NY = MESH.NY; 
	NZ = MESH.NZ;

	StepsPantalla = R1.ProblemNumericalData[8];
	StepsFile = R1.ProblemNumericalData[9];

	ConvergenciaGS = R1.ProblemData[3];
	ConvergenciaGlobal = R1.ProblemData[4];

	EsquemaLargo = R1.ConvectiveScheme1;
	EsquemaCorto = R1.ConvectiveScheme2;

	DIRECTORIO = InputDirectorio;

	//Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	Xcentroide = R1.GeometryData[3];
	Ycentroide = R1.GeometryData[4];

	Xcuadrado = R1.GeometryData[5];
	Ycuadrado = R1.GeometryData[6];

	//Parámetros de computación paralela
	Rank = MPI1.Rank;
	Procesos = MPI1.Procesos;
	Ix = MESH.Ix;
	Fx = MESH.Fx;
	Halo = MPI1.Halo;

	HaloPressure = 1;
	HaloU = 2;
	HaloV = 2;
	HaloW = 2;

	HP = MESH.HP;
	
	//Datos Físicos del Problema
	Rho = R1.ProblemPhysicalData[0];
	Uref = R1.ProblemPhysicalData[1];
	Reynolds = R1.ProblemPhysicalData[2];

	Rayleigh = R1.ProblemPhysicalData[3];
	Cp = R1.ProblemPhysicalData[4];
	Prandtl = R1.ProblemPhysicalData[5];

	gx = R1.ProblemPhysicalData[6];
	gy = R1.ProblemPhysicalData[7];
	gz = R1.ProblemPhysicalData[8];

	//Tleft = R1.ProblemPhysicalData[9]; 
	//Tright = R1.ProblemPhysicalData[10]; 

	//Tbot = R1.ProblemPhysicalData[11]; 
	//Ttop = R1.ProblemPhysicalData[12]; 

	//Variables del Runge - Kutta de 4o order (Capuano et All)
	
	c1 = 0.0;
	c3 = 0.25;
	c2 = (c3 - 1.0)/(4.0*c3 - 3.0);
	c4 = 1.0;

	b1 = 1.0/(12.0*(c3 - 1.0));
	b2 = pow(4.0*c3 - 3.0,2.0)/(12.0*(c3 - 1.0)*(2.0*c3 - 1.0));
	b3 = - 1.0/(12.0*(c3 - 1.0)*(2.0*c3 - 1.0));
	b4 = (4.0*c3 - 3.0)/(12.0*(c3 - 1.0));


	a_21 = (c3 - 1.0)/(4.0*c3 - 3.0);

	a_31 = c3 - ((2.0*c3 - 1.0)*(4.0*c3 - 3.0))/(2.0*(c3 - 1));
	a_32 = ((2.0*c3 - 1.0)*(4.0*c3 - 3.0))/(2.0*(c3 - 1));

	a_41 = - pow(2.0*c3 - 1.0,2.0)/(2.0*(c3 - 1.0)*(4.0*c3 - 3.0));
	a_42 = (6.0*pow(c3,2.0) - 8.0*c3 + 3.0)/(2.0*(c3 - 1.0)*(2.0*c3 - 1.0));
	a_43 = (c3 - 1.0)/((2.0*c3 - 1.0)*(4.0*c3 - 3.0));

/*
	c1 = 0.0;
	c2 = 0.5;
	c3 = 0.5;
	c4 = 1.0;

	b1 = 1.0/6.0;
	b2 = 1.0/3.0;
	b3 = 1.0/3.0;
	b4 = 1.0/6.0;


	a_21 = 0.5;

	a_31 = 0.0;
	a_32 = 0.5;

	a_41 = 0.0;
	a_42 = 0.0;
	a_43 = 1.0;
*/

	//Cálculos extra para cada problema
	if(Problema == 1){ //Problema Driven Cavity

		mu = (Rho*Uref*Xdominio)/Reynolds;

	}

	else if(Problema == 2){ //Problema DIfferentially Heated

		//Cálculo de To
		if(abs(0.50*(Tleft + Tright)) > 0.0){
			To = abs(0.50*(Tleft + Tright));
			Difference = abs(Tleft - Tright);
		}
		else{
			To = abs(0.50*(Tbot + Ttop));
			Difference = abs(Tbot - Ttop);
		}

		Beta = 1.0/To;
		Producto = (pow(Rho,2)*abs(gy)*pow(Xdominio,3)*Beta*Difference*Prandtl)/Rayleigh;

		mu = sqrt(Producto);

		K = (Cp*mu)/Prandtl;

	}

}

//Alojamiento de memoria para las matrices necesarias
void Solver::AllocateMatrix(Memory M1){

	if(Rank == 0){

		//Matrices globales de propiedades
		PGlobal = M1.AllocateDouble(NX + 2*HaloPressure, NY + 2*HaloPressure, NZ + 2*HaloPressure, 1); //Presión P Global
		UGlobal = M1.AllocateDouble(NX + 1 + 2*HaloU, NY + 2*HaloU, NZ + 2*HaloU, 1); //Velocidad U Global
		VGlobal = M1.AllocateDouble(NX + 2*HaloV, NY + 1 + 2*HaloV, NZ + 2*HaloV, 1); //Velocidad V Global
		WGlobal = M1.AllocateDouble(NX + 2*HaloW, NY + 2*HaloW, NZ + 1 + 2*HaloW, 1); //Velocidad W Global

	}
	
}





//Actualización del Halo de Velocidades en cada Step
void Solver::Get_HaloVelocities(){
int i, j, k;

	//Velocidad U
	if(Rank == 0){

		//Halo Izquierda
		for(i = - HaloU; i < 0; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LU(i,j,k,0)] = Uleft[ULEFT(0,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HP; j++){
					ULPRES[LU(i,j,k,0)] = Utop[UTOP(0,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HP; j < 0; j++){
					ULPRES[LU(i,j,k,0)] = Ubot[UBOT(0,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LU(i,j,k,0)] = Uhere[UHERE(0,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HP; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LU(i,j,k,0)] = Uthere[UTHERE(0,j,0)];
				}
			}

		}

	}
	else if(Rank == Procesos - 1){

		//Halo Derecha
		for(i = NX + 1; i < NX + HaloU + 1; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LU(i,j,k,0)] = Uright[URIGHT(NX,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HP; j++){
					ULPRES[LU(i,j,k,0)] = Utop[UTOP(NX,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HP; j < 0; j++){
					ULPRES[LU(i,j,k,0)] = Ubot[UBOT(NX,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LU(i,j,k,0)] = Uhere[UHERE(NX,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HP; k++){
				for(j = 0; j < NY; j++){
					ULPRES[LU(i,j,k,0)] = Uthere[UTHERE(NX,j,0)];
				}
			}

		}

	}

	//Resto de Halos
	for(i = Ix; i < Fx + 1; i++){

		//Halo Superior
		for(k = 0; k < NZ; k++){
			for(j = NY; j < NY + HaloU; j++){
				ULPRES[LU(i,j,k,0)] = Utop[UTOP(i,NY,k)];
			}
		}

		//Halo Inferior
		for(k = 0; k < NZ; k++){
			for(j = - HaloU; j < 0; j++){
				ULPRES[LU(i,j,k,0)] = Ubot[UBOT(i,0,k)];
			}
		}

		//Halo Here
		for(k = - HaloU; k < 0; k++){
			for(j = 0; j < NY; j++){
				ULPRES[LU(i,j,k,0)] = ULPRES[LU(i,j,NZ + k,0)];
			}
		}

		//Halo There
		for(k = NZ; k < NZ + HaloU; k++){
			for(j = 0; j < NY; j++){
				ULPRES[LU(i,j,k,0)] = ULPRES[LU(i,j,k - NZ,0)];
			}
		}

	}	

	//Velocidad V
	if(Rank == 0){

		//Halo Izquierda
		for(i = - HaloV; i < 0; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 1; j < NY; j++){
					VLPRES[LV(i,j,k,0)] = Vleft[VLEFT(0,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HaloV + 1; j++){
					VLPRES[LV(i,j,k,0)] = Vtop[VTOP(0,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HaloV; j <= 0; j++){
					VLPRES[LV(i,j,k,0)] = Vbot[VBOT(0,0,k)];
				}
			}

			//Halo Here
			for(k = - HP; k < 0; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LV(i,j,k,0)] = Vhere[VHERE(0,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HaloV; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LV(i,j,k,0)] = Vthere[VTHERE(0,j,0)];
				}
			}

		}

	}
	else if(Rank == Procesos - 1){

		//Halo Derecha
		for(i = NX; i < NX + HaloV; i++){

			//Halo Channel
			for(k = 0; k < NZ; k++){
				for(j = 1; j < NY; j++){
					VLPRES[LV(i,j,k,0)] = Vright[VRIGHT(NX,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ; k++){
				for(j = NY; j < NY + HaloV; j++){
					VLPRES[LV(i,j,k,0)] = Vtop[VTOP(NX,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ; k++){
				for(j = - HaloV; j <= 0; j++){
					VLPRES[LV(i,j,k,0)] = Vbot[VBOT(NX,0,k)];
				}
			}

			//Halo Here
			for(k = - HaloV; k < 0; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LV(i,j,k,0)] = Vhere[VHERE(NX,j,0)];
				}
			}

			//Halo There
			for(k = NZ; k < NZ + HaloV; k++){
				for(j = 0; j < NY + 1; j++){
					VLPRES[LV(i,j,k,0)] = Vthere[VTHERE(NX,j,0)];
				}
			}

		}

	}

	//Resto de Halos
	for(i = Ix - 1; i < Fx + 1; i++){

		//Halo Superior
		for(k = 0; k < NZ; k++){
			for(j = NY + 1; j < NY + HaloV + 1; j++){
				VLPRES[LV(i,j,k,0)] = Vtop[VTOP(i,NY,k)];
			}
		}

		//Halo Inferior
		for(k = 0; k < NZ; k++){
			for(j = - HaloV; j < 0; j++){
				VLPRES[LV(i,j,k,0)] = Vbot[VBOT(i,0,k)];
			}
		}

		//Halo Here
		for(k = - HaloV; k < 0; k++){
			for(j = 0; j < NY + 1; j++){
				VLPRES[LV(i,j,k,0)] = VLPRES[LV(i,j,NZ + k,0)];
			}
		}

		//Halo There
		for(k = NZ; k < NZ + HaloV; k++){
			for(j = 0; j < NY + 1; j++){
				VLPRES[LV(i,j,k,0)] = VLPRES[LV(i,j,k - NZ,0)];
			}
		}

	}


	//Velocidad W
	if(Rank == 0){

		//Halo Izquierda
		for(i = - HaloW; i < 0; i++){

			//Halo Channel
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LW(i,j,k,0)] = Wleft[WLEFT(0,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ + 1; k++){
				for(j = NY; j < NY + HaloW; j++){
					WLPRES[LW(i,j,k,0)] = Wtop[WTOP(0,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ + 1; k++){
				for(j = - HaloW; j < 0; j++){
					WLPRES[LW(i,j,k,0)] = Wbot[WBOT(0,0,k)];
				}
			}

			//Halo Here
			for(k = - HaloW; k < 0; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LW(i,j,k,0)] = Where[WHERE(0,j,0)];
				}
			}

			//Halo There
			for(k = NZ + 1; k < NZ + HaloW + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LW(i,j,k,0)] = Wthere[WTHERE(0,j,0)];
				}
			}

		}

	}
	else if(Rank == Procesos - 1){

		//Halo Derecha
		for(i = NX; i < NX + HaloW; i++){

			//Halo Channel
			for(k = 0; k < NZ + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LW(i,j,k,0)] = Wright[WRIGHT(NX,j,k)];
				}
			}

			//Halo Superior
			for(k = 0; k < NZ + 1; k++){
				for(j = NY; j < NY + HaloW; j++){
					WLPRES[LW(i,j,k,0)] = Wtop[WTOP(NX,0,k)];
				}
			}

			//Halo Inferior
			for(k = 0; k < NZ + 1; k++){
				for(j = - HaloW; j < 0; j++){
					WLPRES[LW(i,j,k,0)] = Wbot[WBOT(NX,0,k)];
				}
			}

			//Halo Here
			for(k = - HaloW; k < 0; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LW(i,j,k,0)] = Where[WHERE(NX,j,0)];
				}
			}

			//Halo There
			for(k = NZ + 1; k < NZ + HaloW + 1; k++){
				for(j = 0; j < NY; j++){
					WLPRES[LW(i,j,k,0)] = Wthere[WTHERE(NX,j,0)];
				}
			}

		}

	}

	//Resto de Halos
	for(i = Ix - 1; i < Fx + 1; i++){

		//Halo Superior
		for(k = 0; k < NZ + 1; k++){
			for(j = NY; j < NY + HaloW; j++){
				WLPRES[LW(i,j,k,0)] = Wtop[WTOP(i,NY,k)];
			}
		}

		//Halo Inferior
		for(k = 0; k < NZ + 1; k++){
			for(j = - HaloW; j < 0; j++){
				WLPRES[LW(i,j,k,0)] = Wbot[WBOT(i,0,k)];
			}
		}

		//Halo Here
		for(k = - HaloW; k < 0; k++){
			for(j = 0; j < NY; j++){
				WLPRES[LW(i,j,k,0)] = WLPRES[LW(i,j,NZ + k,0)];
			}
		}

		//Halo There
		for(k = NZ + 1; k < NZ + HaloW + 1; k++){
			for(j = 0; j < NY; j++){
				WLPRES[LW(i,j,k,0)] = WLPRES[LW(i,j,k - NZ,0)];
			}
		}

	}
	//Problema Differentially Heated
	if(Problema == 2){

		if(Rank == 0){

			//Parte Izquierda
			for(i = - HP; i < 0; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){
						TLPRES[LP(i,j,k,0)] = TLEFT[PLEFT(0,j,k)];
					}
				}
			}

		}
		else if(Rank == Procesos - 1){

			//Parte Derecha
			for(i = NX; i < NX + HP; i++){
				for(j = 0; j < NY; j++){
					for(k = 0; k < NZ; k++){
						TLPRES[LP(i,j,k,0)] = TRIGHT[PRIGHT(0,j,k)];
					}
				}
			}

		}

		//Resto de Halos
	
		for(i = Ix; i < Fx; i++){

			//Partes Superior e Inferior
			for(k = 0; k < NZ; k++){

				//Parte Superior
				for(j = NY; j < NY + HP; j++){
					TLPRES[LP(i,j,k,0)] = TTOP[PTOP(i,NY,k)];
				}
				
				//Parte Inferior
				for(j = - HP; j < 0; j++){
					TLPRES[LP(i,j,k,0)] = TBOT[PBOT(i,0,k)];
				}

			}

			//Partes Here y There
			for(j = 0; j < NY; j++){

				//Parte Here
				for(k = - HP; k < 0; k++){
					TLPRES[LP(i,j,k,0)] = There[PHERE(i,j,0)];
				}

				//Parte There
				for(k = NZ; k < NZ + HP; k++){
					TLPRES[LP(i,j,k,0)] = Tthere[PTHERE(i,j,NZ)];
				}

			}

		}

	}

}





void Solver::Get_PressureCoefficients(Mesher MESH){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Parte Central
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){		
					aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
					ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

					as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
					an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

					ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
					at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
				}
			}		
		}

		//Parte Inferior
		j = 0;

		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = 0.0;
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte Superior
		j = NY - 1;

		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){	
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = 0.0;

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte Here
		k = 0;

		for(i = Ix; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte There
		k = NZ - 1;

		for(i = Ix; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);
			}		
		}

		//Esquinas
		for(i = Ix; i < Fx; i++){

			//Esquina Abajo Here
			j = 0;
			k = 0;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
		
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);


			//Esquina Arriba There
			j = NY - 1; 
			k = NZ - 1;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);
		}

	}
	else if(Rank == 0){

		//Parte Central
		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){		
					aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
					ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

					as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
					an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

					ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
					at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
				}
			}		
		}

		//Parte Inferior
		j = 0;

		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = 0.0;
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];

			}		
		}

		//Parte Superior
		j = NY - 1;

		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){	
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = 0.0;

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte Here
		k = 0;

		for(i = Ix + 1; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte There
		k = NZ - 1;

		for(i = Ix + 1; i < Fx; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);
			}		
		}

		//Esquinas
		for(i = Ix + 1; i < Fx; i++){

			//Esquina Abajo Here
			j = 0;
			k = 0;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
		
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);


			//Esquina Arriba There
			j = NY - 1; 
			k = NZ - 1;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);

		}

		//Parte Izquierda
		i = 0;

		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = 0.0;
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}
		}

		//Parte Inferior
		j = 0;

		for(k = 1; k < NZ - 1; k++){		
			aw[LA(i,j,k,0)] = 0.0;
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
		}		
		
		//Parte Superior
		j = NY - 1;

		for(k = 1; k < NZ - 1; k++){	
			aw[LA(i,j,k,0)] = 0.0;
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];	
		}

		//Parte Here
		k = 0;

		for(j = 1; j < NY - 1; j++){		
			aw[LA(i,j,k,0)] = 0.0;
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];	
		}

		//Parte There
		k = NZ - 1;

		for(j = 1; j < NY - 1; j++){		
			aw[LA(i,j,k,0)] = 0.0;
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);
		}

		//Esquina Abajo Here
		j = 0;
		k = 0;

		aw[LA(i,j,k,0)] = 0.0;
		ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

		as[LA(i,j,k,0)] = 0.0;
		an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
		
		//Esquina Arriba Here
		j = NY - 1;
		k = 0;

		aw[LA(i,j,k,0)] = 0.0;
		ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

		as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
		an[LA(i,j,k,0)] = 0.0;

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];

		//Esquina Abajo There
		j = 0;
		k = NZ - 1;

		aw[LA(i,j,k,0)] = 0.0;
		ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

		as[LA(i,j,k,0)] = 0.0;
		an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);


		//Esquina Arriba There
		j = NY - 1; 
		k = NZ - 1;

		aw[LA(i,j,k,0)] = 0.0;
		ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

		as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
		an[LA(i,j,k,0)] = 0.0;

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);


	}
	else if(Rank == Procesos - 1){

		//Parte Central
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){		
					aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
					ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

					as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
					an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

					ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
					at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
				}
			}		
		}

		//Parte Inferior
		j = 0;

		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = 0.0;
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte Superior
		j = NY - 1;

		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){	
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = 0.0;

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte Here
		k = 0;

		for(i = Ix; i < Fx - 1; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte There
		k = NZ - 1;

		for(i = Ix; i < Fx - 1; i++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);
			}		
		}

		//Esquinas
		for(i = Ix; i < Fx - 1; i++){

			//Esquina Abajo Here
			j = 0;
			k = 0;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
		
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);


			//Esquina Arriba There
			j = NY - 1; 
			k = NZ - 1;

			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,1)]/MESH.DeltasMU[GU(i+1,j,k,0)];

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);

		}


		//Parte Derecha
		i = NX - 1;

		//Parte Central
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY - 1; j++){		
				aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
				ae[LA(i,j,k,0)] = 0.0;

				as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
				an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

				ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
				at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
			}		
		}

		//Parte Inferior
		j = 0;

		for(k = 1; k < NZ - 1; k++){		
			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = 0.0;

			as[LA(i,j,k,0)] = 0.0;
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];		
		}

		//Parte Superior
		j = NY - 1;

		for(k = 1; k < NZ - 1; k++){	
			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = 0.0;

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = 0.0;

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];
		
		}

		//Parte Here
		k = 0;

		for(j = 1; j < NY - 1; j++){		
			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = 0.0;

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];		
		}

		//Parte There
		k = NZ - 1;

		for(j = 1; j < NY - 1; j++){		
			aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
			ae[LA(i,j,k,0)] = 0.0;

			as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
			an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

			ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
			at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);
		}

		//Esquinas
		
		//Esquina Abajo Here
		j = 0;
		k = 0;

		aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
		ae[LA(i,j,k,0)] = 0.0;

		as[LA(i,j,k,0)] = 0.0;
		an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];	
		
		//Esquina Arriba Here
		j = NY - 1;
		k = 0;

		aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
		ae[LA(i,j,k,0)] = 0.0;

		as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
		an[LA(i,j,k,0)] = 0.0;

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/(2.0*MESH.DeltasMW[GW(i,j,k,2)]);
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/MESH.DeltasMW[GW(i,j,k+1,2)];	

		//Esquina Abajo There
		j = 0;
		k = NZ - 1;

		aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
		ae[LA(i,j,k,0)] = 0.0;

		as[LA(i,j,k,0)] = 0.0;
		an[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,3)]/MESH.DeltasMV[GV(i,j+1,k,1)];

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);


		//Esquina Arriba There
		j = NY - 1; 
		k = NZ - 1;

		aw[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,0)]/MESH.DeltasMU[GU(i,j,k,0)];
		ae[LA(i,j,k,0)] = 0.0;

		as[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,2)]/MESH.DeltasMV[GV(i,j,k,1)];
		an[LA(i,j,k,0)] = 0.0;

		ah[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,4)]/MESH.DeltasMW[GW(i,j,k,2)];
		at[LA(i,j,k,0)] = MESH.SupMP[GP(i,j,k,5)]/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)]);

	}

	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
				ap[LA(i,j,k,0)] = aw[LA(i,j,k,0)] + ae[LA(i,j,k,0)] + as[LA(i,j,k,0)] + an[LA(i,j,k,0)] + ah[LA(i,j,k,0)] + at[LA(i,j,k,0)];
			}
		}
	}

}

//Cálculo del término difusivo de la velocidad U
void Solver::Get_DiffusiveU(Mesher MESH, double *UFIELD){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx + 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		
		for(i = Ix; i < Fx + 1; i++){

			//Partes Inferior y Superior
			for(k = 1; k < NZ - 1; k++){
				
				//Parte Inferior
				j = 0;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}
			
			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Flas de las Esquinas

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			
		}
	
	}
	else if(Rank == 0){

		//Centro
		for(i = Ix + 1; i < Fx + 1; i++){
			for(k = 1; k < NZ-1; k++){
				for(j = 1; j < NY-1; j++){
					DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix + 1; i < Fx + 1; i++){

			//Partes Inferior y Superior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Flas de las Esquinas

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

	}
	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		for(i = Ix; i < Fx; i++){

			//Partes Inferior y Superior
			for(k = 1; k < NZ - 1; k++){

				//Parte Inferior
				j = 0;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

				//Parte Superior
				j = NY - 1;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Partes Here y There
			for(j = 1; j < NY - 1; j++){

				//Parte Here
				k = 0;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

				//Parte There
				k = NZ - 1;

				DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

			//Filas de las Esquinas

			//Fila Abajo Here
			j = 0;
			k = 0;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

			//Fila Abajo There
			j = 0;
			k = NZ - 1;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(UFIELD[LU(i,j+1,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - Ubot[UBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			//Fila Arriba Here
			j = NY - 1;
			k = 0;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,k+1,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

			//Fila Arriba There
			j = NY - 1;
			k = NZ - 1;

			DiffusiveU[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[GU(i,j,k,0)]))*(
											+ MESH.SupMU[GU(i,j,k,1)]*(UFIELD[LU(i+1,j,k,0)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,0)]
											- MESH.SupMU[GU(i,j,k,0)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i-1,j,k,0)])/MESH.DeltasMP[GP(i-1,j,k,0)]
											+ MESH.SupMU[GU(i,j,k,3)]*(Utop[UTOP(i,j,k)] - UFIELD[LU(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMU[GU(i,j,k,2)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMU[GU(i,j,k,5)]*(UFIELD[LU(i,j,0,0)] - UFIELD[LU(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMU[GU(i,j,k,4)]*(UFIELD[LU(i,j,k,0)] - UFIELD[LU(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

	}	
		
}

//Cálculo del término difusivo de la velocidad V
void Solver::Get_DiffusiveV(Mesher MESH, double *VFIELD){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}
		}

		//Partes Here y There
		for(i = Ix; i < Fx; i++){

			for(j = 1; j < NY; j++){

				//Parte Here
				k = 0;

				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*( 
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

				//Parte There
				k = NZ - 1;

				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,0,0)] - VFIELD[LV(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}
		}

	}
	else if(Rank == 0){
	
		for(i = Ix + 1; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY; j++){

				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY; j++){
				
				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,0,0)] - VFIELD[LV(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

		}

		
		

		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - Vleft[VLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}
		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){

			DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - Vleft[VLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

		}
		
		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){
		
			DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - Vleft[VLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,0,0)] - VFIELD[LV(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
		}
	
	}
	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;

			for(j = 1; j < NY; j++){

				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);
			}

			//Parte There
			k = NZ - 1;

			for(j = 1; j < NY; j++){

				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(VFIELD[LV(i+1,j,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,0,0)] - VFIELD[LV(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

			}

		}


		//Parte Derecha
		i = NX - 1;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(Vright[VRIGHT(i,j,k)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);
			}
		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){

			DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(Vright[VRIGHT(i,j,k)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,k+1,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMW[GW(i,j,k+1,2)]
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j,NZ-1,0)])/(2.0*MESH.DeltasMW[GW(i,j,k,2)])
											);

		}

		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){
			
			DiffusiveV[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[GV(i,j,k,0)]))*(
											+ MESH.SupMV[GV(i,j,k,1)]*(Vright[VRIGHT(i,j,k)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMV[GV(i,j,k,0)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMV[GV(i,j,k,3)]*(VFIELD[LV(i,j+1,k,0)] - VFIELD[LV(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,1)]
											- MESH.SupMV[GV(i,j,k,2)]*(VFIELD[LV(i,j,k,0)] - VFIELD[LV(i,j-1,k,0)])/MESH.DeltasMP[GP(i,j-1,k,1)]
											+ MESH.SupMV[GV(i,j,k,5)]*(VFIELD[LV(i,j,0,0)] - VFIELD[LV(i,j,k,0)])/(2.0*MESH.DeltasMW[GW(i,j,k+1,2)])
											- MESH.SupMV[GV(i,j,k,4)]*(VFIELD[LV(i,j,k,0)] - VLPRES[LV(i,j,k-1,0)])/MESH.DeltasMW[GW(i,j,k,2)]
											);

		}

	}

}

//Cálculo del término difusivo de la velocidad W
void Solver::Get_DiffusiveW(Mesher MESH, double *WFIELD){
int i, j, k;

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro
		for(i = Ix; i < Fx; i++){
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}

		}
		
	}
	else if(Rank == 0){

		//Centro
		for(i = Ix + 1; i < Fx; i++){
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}

		}

		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){
				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - Wleft[WLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ; k++){

			DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - Wleft[WLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){

			DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - Wleft[WLEFT(i,j,k)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}

	}
	else if(Rank == Procesos - 1){

		//Centro
		for(i = Ix; i < Fx - 1; i++){
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){
					DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);
				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(WFIELD[LW(i+1,j,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

			}

		}

		//Parte Derecha
		i = NX - 1;

		//Centro
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){

				DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*( 0.0
											+ MESH.SupMW[GW(i,j,k,1)]*(Wright[WRIGHT(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;

		for(k = 1; k < NZ; k++){

			DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(Wright[WRIGHT(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(WFIELD[LW(i,j+1,k,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - Wbot[WBOT(i,j,k)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){
			
			DiffusiveW[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[GW(i,j,k,0)]))*(
											+ MESH.SupMW[GW(i,j,k,1)]*(Wright[WRIGHT(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMU[GU(i+1,j,k,0)]
											- MESH.SupMW[GW(i,j,k,0)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i-1,j,k,0)])/MESH.DeltasMU[GU(i,j,k,0)]
											+ MESH.SupMW[GW(i,j,k,3)]*(Wtop[WTOP(i,j,k)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMV[GV(i,j+1,k,1)]
											- MESH.SupMW[GW(i,j,k,2)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j-1,k,0)])/MESH.DeltasMV[GV(i,j,k,1)]
											+ MESH.SupMW[GW(i,j,k,5)]*(WFIELD[LW(i,j,k+1,0)] - WFIELD[LW(i,j,k,0)])/MESH.DeltasMP[GP(i,j,k,2)]
											- MESH.SupMW[GW(i,j,k,4)]*(WFIELD[LW(i,j,k,0)] - WFIELD[LW(i,j,k-1,0)])/MESH.DeltasMP[GP(i,j,k-1,2)]
											);

		}
	
	}

}

//Cálculo del término convectivo de la velocidad U
void Solver::Get_ConvectiveU(Mesher MESH, double *UFIELD, double *VFIELD, double *WFIELD){
int i, j, k;
double uW, uE, uS, uN, uH, uT, vS, vN, wH, wT;
double uW_pred, uE_pred, uS_pred, uN_pred, uH_pred, uT_pred, vS_pred, vN_pred, wH_pred, wT_pred;

	//ESQUEMA CONVECTIVO:
	//Corto (Interpolación)
	//(CoordObjetivo, Coord1, Valor1, Coord2, Valor2)


	//Largo
	//(CoordObjetivo, Velocidad, Coord1, Valor1, Coord2, Valor2, Coord3, Valor3, Coord4, Valor4, EsquemaLargo)

	if(Rank != 0 && Rank != Procesos - 1){

		//Centro 
		for(i = Ix; i < Fx + 1; i++){
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], ULPRES[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

					uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
					uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

					uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
					uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

					vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

					uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

					uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
					uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

					uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
					uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
					wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

					ConvectiveU[LU(i,j,k,0)] = (1.0/MESH.VolMU[GU(i,j,k,0)])*(
											  - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											  + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											  - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											  + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											  - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											  + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											  );
				}	
			}

			//Partes Inferior y Superior	

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Ubot[UBOT(i,j,k)];
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = Ubot[UBOT(i,j,k)];
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Utop[UTOP(i,j,k)];

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = Utop[UTOP(i,j,k)];

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}


			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j <  NY - 1; j++){

				
				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Uhere[UHERE(i,j,k)];
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = Uhere[UHERE(i,j,k)];
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Uthere[UTHERE(i,j,k)];

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = Uthere[UTHERE(i,j,k)];

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Esquina Abajo Here
			j = 0;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

		}		

	}
	else if(Rank == 0){

		 
		for(i = Ix + 1; i < Fx + 1; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

					uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
					uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

					uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
					uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

					vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

					uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

					uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
					uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

					uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
					uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
					wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

					ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Ubot[UBOT(i,j,k)];
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = Ubot[UBOT(i,j,k)];
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Utop[UTOP(i,j,k)];

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = Utop[VTOP(i,j,k)];

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}


			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j <  NY - 1; j++){	

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Uhere[UHERE(i,j,k)];
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = Uhere[UHERE(i,j,k)];
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Uthere[UTHERE(i,j,k)];

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = Uthere[UTHERE(i,j,k)];

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Esquina Abajo Here
			j = 0;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );


			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );


		}	

	}
	else if(Rank == Procesos - 1){

		 
		for(i = Ix; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY - 1; j++){

					uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

					uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
					uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

					uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
					uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

					vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

					uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

					uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
					uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

					uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
					uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
					wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

					ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

				}
			}


			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Ubot[UBOT(i,j,k)];
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = Ubot[UBOT(i,j,k)];
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;

			for(k = 1; k < NZ - 1; k++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Utop[UTOP(i,j,k)];

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = Utop[UTOP(i,j,k)];

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Partes Here y There

			//Parte Here
			k = 0;

			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Uhere[UHERE(i,j,k)];
				uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = Uhere[UHERE(i,j,k)];
				uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Parte There
			k = NZ - 1;

			for(j = 1; j <  NY - 1; j++){

				uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

				uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

				uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uT_pred = Uthere[UTHERE(i,j,k)];

				vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

				uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

				uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uT = Uthere[UTHERE(i,j,k)];

				vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
				wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			}

			//Esquina Abajo Here
			j = 0;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Ubot[UBOT(i,j,k)];
			uN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)]);

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = Ubot[UBOT(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], MESH.MU[GU(i,j+2,k,1)], UFIELD[LU(i,j+2,k,0)], EsquemaLargo);

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = 0.50*(Vbot[VBOT(i - 1,j,k)] + Vbot[VBOT(i,j,k)]);
			vN = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uN_pred, MESH.MV[GV(i-2,j+1,k,0)], VFIELD[LV(i-2,j+1,k,0)], MESH.MV[GV(i-1,j+1,k,0)], VFIELD[LV(i-1,j+1,k,0)], MESH.MV[GV(i,j+1,k,0)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i+1,j+1,k,0)], VFIELD[LV(i+1,j+1,k,0)], EsquemaLargo);
					
			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba Here
			j = NY - 1;
			k = 0;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Uhere[UHERE(i,j,k)];
			uT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)]);

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = Uhere[UHERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], MESH.MU[GU(i,j,k+2,2)], UFIELD[LU(i,j,k+2,0)], EsquemaLargo);

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = 0.50*(Where[WHERE(i - 1,j,k)] + Where[WHERE(i,j,k)]);
			wT = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uT_pred, MESH.MW[GW(i-2,j,k+1,0)], WFIELD[LW(i-2,j,k+1,0)], MESH.MW[GW(i-1,j,k+1,0)], WFIELD[LW(i-1,j,k+1,0)], MESH.MW[GW(i,j,k+1,0)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i+1,j,k+1,0)], WFIELD[LW(i+1,j,k+1,0)], EsquemaLargo);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;

			uW_pred = Interpolacion(MESH.MP[GP(i-1,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)]);
			uE_pred = Interpolacion(MESH.MP[GP(i,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)]);

			uS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
			uN_pred = Utop[UTOP(i,j,k)];

			uH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
			uT_pred = Uthere[UTHERE(i,j,k)];

			vS_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
			vN_pred = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
			wT_pred = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			uW = ConvectiveScheme(MESH.MP[GP(i-1,j,k,0)], uW_pred, MESH.MU[GU(i-2,j,k,0)], UFIELD[LU(i-2,j,k,0)], MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[GP(i,j,k,0)], uE_pred, MESH.MU[GU(i-1,j,k,0)], UFIELD[LU(i-1,j,k,0)], MESH.MU[GU(i,j,k,0)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i+1,j,k,0)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+2,j,k,0)], UFIELD[LU(i+2,j,k,0)], EsquemaLargo);

			uS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
			uN = Utop[UTOP(i,j,k)];

			uH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
			uT = Uthere[UTHERE(i,j,k)];

			vS = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uS_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
			vN = 0.50*(Vtop[VTOP(i - 1,j,k)] + Vtop[VTOP(i,j,k)]);

			wH = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uH_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
			wT = 0.50*(Wthere[WTHERE(i - 1,j,k)] + Wthere[WTHERE(i,j,k)]);

			ConvectiveU[LU(i,j,k,0)] = (1/MESH.VolMU[GU(i,j,k,0)])*(
											 - MESH.SupMU[GU(i,j,k,0)]*uW*uW
											 + MESH.SupMU[GU(i,j,k,1)]*uE*uE
											 - MESH.SupMU[GU(i,j,k,2)]*uS*vS
											 + MESH.SupMU[GU(i,j,k,3)]*uN*vN
											 - MESH.SupMU[GU(i,j,k,4)]*uH*wH
											 + MESH.SupMU[GU(i,j,k,5)]*uT*wT
											 );

		}	
	}
}

//Cálculo del término convectivo de la velocidad V
void Solver::Get_ConvectiveV(Mesher MESH, double *UFIELD, double *VFIELD, double *WFIELD){
int i, j, k;
double vW, vE, vS, vN, vH, vT, uW, uE, wH, wT;
double vW_pred, vE_pred, vS_pred, vN_pred, vH_pred, vT_pred, uW_pred, uE_pred, wH_pred, wT_pred;

	//ESQUEMA CONVECTIVO:
	//Corto (Interpolación)
	//(CoordObjetivo, Coord1, Valor1, Coord2, Valor2)


	//Largo
	//(CoordObjetivo, Velocidad, Coord1, Valor1, Coord2, Valor2, Coord3, Valor3, Coord4, Valor4, EsquemaLargo)

	if(Rank != 0 && Rank != Procesos - 1){

		for(i = Ix; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
					vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

					vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
					vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

					uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

					wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

					vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
					vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

					vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
					vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				}
			}

			//Partes Here y There
			//Parte Here
			k = 0;

			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

			
				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

				
				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				
			}

		}
		
	}	
	else if(Rank == 0){

		
		for(i = Ix + 1; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
					vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

					vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
					vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

					uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

					wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

					vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
					vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

					vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
					vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;
			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);


				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY; j++){
				
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}	

		}

		//Parte Izquierda
		i = 0;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){

				vW_pred = Vleft[VLEFT(i,j,k)];
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

				vW = Vleft[VLEFT(i,j,k)];
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
			}
		}
		
		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){

				

				vW_pred = Vleft[VLEFT(i,j,k)];
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

				vW = Vleft[VLEFT(i,j,k)];
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

				
				

				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
		}
		
		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){
					
				vW_pred = Vleft[VLEFT(i,j,k)];
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = Vleft[VLEFT(i,j,k)];
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = 0.50*(Uleft[ULEFT(i,j - 1,k)] + Uleft[ULEFT(i,j,k)]);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);


				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

		}

	}
	else if(Rank == Procesos - 1){
		
		
		for(i = Ix; i < Fx - 1; i++){

			//Centro
			for(k = 1; k < NZ - 1; k++){
				for(j = 1; j < NY; j++){
					vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
					vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

					vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
					vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

					uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

					wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

					vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
					vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

					vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
					vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				}
			}

			//Partes Here y There

			//Parte Here
			k = 0;

			for(j = 1; j < NY; j++){
				
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,0)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
			}

			//Parte There
			k = NZ - 1;

			for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], MESH.MV[GV(i+2,j,k,0)], VFIELD[LV(i+2,j,k,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vE_pred, MESH.MU[GU(i+1,j-2,k,1)], UFIELD[LU(i+1,j-2,k,0)], MESH.MU[GU(i+1,j-1,k,1)], UFIELD[LU(i+1,j-1,k,0)], MESH.MU[GU(i+1,j,k,1)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j+1,k,1)], UFIELD[LU(i+1,j+1,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,0)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );

			}
		}

		//Parte Derecha
		i = NX - 1;

		//Centro
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Vright[VRIGHT(i,j,k)];

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = Vright[VRIGHT(i,j,k)];

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
			}
		}

		//Partes Here y There

		//Parte Here
		k = 0;
		for(j = 1; j < NY; j++){
				
				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Vright[VRIGHT(i,j,k)];

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Vhere[VHERE(i,j,k)];
				vT_pred = Interpolacion(MESH.MW[GW(i,j,k+1,2)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,2)]);

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH_pred = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = Vright[VRIGHT(i,j,k)];

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = Vhere[VHERE(i,j,k)];
				vT = ConvectiveScheme(MESH.MW[GW(i,j,k+1,2)], wT_pred, MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], MESH.MV[GV(i,j,k+2,2)], VFIELD[LV(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH = 0.50*(Where[WHERE(i,j - 1,k)] + Where[WHERE(i,j,k)]);
				wT = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vT_pred, MESH.MW[GW(i,j-2,k+1,1)], WFIELD[LW(i,j-2,k+1,0)], MESH.MW[GW(i,j-1,k+1,1)], WFIELD[LW(i,j-1,k+1,0)], MESH.MW[GW(i,j,k+1,1)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j+1,k+1,1)], WFIELD[LW(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*(  
											 - MESH.SupMV[GV(i,j,k,0)]*uW*vW
											 + MESH.SupMV[GV(i,j,k,1)]*uE*vE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
		}


		//Parte There
		k = NZ - 1;
		for(j = 1; j < NY; j++){

				vW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)]);
				vE_pred = Vright[VRIGHT(i,j,k)];

				vS_pred = Interpolacion(MESH.MP[GP(i,j-1,k,1)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[GP(i,j,k,1)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)]);

				vH_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,2)]);
				vT_pred = Vthere[VTHERE(i,j,k)];

				uW_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wT_pred = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				vW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MV[GV(i-2,j,k,0)], VFIELD[LV(i-2,j,k,0)], MESH.MV[GV(i-1,j,k,0)], VFIELD[LV(i-1,j,k,0)], MESH.MV[GV(i,j,k,0)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i+1,j,k,0)], VFIELD[LV(i+1,j,k,0)], EsquemaLargo);
				vE = Vright[VRIGHT(i,j,k)];

				vS = ConvectiveScheme(MESH.MP[GP(i,j-1,k,1)], vS_pred, MESH.MV[GV(i,j-2,k,1)], VFIELD[LV(i,j-2,k,0)], MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[GP(i,j,k,1)], vN_pred, MESH.MV[GV(i,j-1,k,1)], VFIELD[LV(i,j-1,k,0)], MESH.MV[GV(i,j,k,1)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j+1,k,1)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+2,k,1)], VFIELD[LV(i,j+2,k,0)], EsquemaLargo);

				vH = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wH_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,2)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vT = Vthere[VTHERE(i,j,k)];

				uW = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vW_pred, MESH.MU[GU(i,j-2,k,1)], UFIELD[LU(i,j-2,k,0)], MESH.MU[GU(i,j-1,k,1)], UFIELD[LU(i,j-1,k,0)], MESH.MU[GU(i,j,k,1)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j+1,k,1)], UFIELD[LU(i,j+1,k,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j - 1,k)] + Uright[URIGHT(i,j,k)]);

				wH = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vH_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wT = 0.50*(Wthere[WTHERE(i,j - 1,k)] + Wthere[WTHERE(i,j,k)]);

				ConvectiveV[LV(i,j,k,0)] = (1.0/MESH.VolMV[GV(i,j,k,0)])*( 
											 - MESH.SupMV[GV(i,j,k,0)]*vW*uW
											 + MESH.SupMV[GV(i,j,k,1)]*vE*uE
											 - MESH.SupMV[GV(i,j,k,2)]*vS*vS
											 + MESH.SupMV[GV(i,j,k,3)]*vN*vN
											 - MESH.SupMV[GV(i,j,k,4)]*vH*wH
											 + MESH.SupMV[GV(i,j,k,5)]*vT*wT
											 );
				
		}

	}

}

//Cálculo del término convectivo de la velocidad W
void Solver::Get_ConvectiveW(Mesher MESH, double *UFIELD, double *VFIELD, double *WFIELD){
int i, j, k;
double wW, wE, wS, wN, wH, wT, uW, uE, vS, vN;
double wW_pred, wE_pred, wS_pred, wN_pred, wH_pred, wT_pred, uW_pred, uE_pred, vS_pred, vN_pred;

	//ESQUEMA CONVECTIVO:
	//Corto (Interpolación)
	//(CoordObjetivo, Coord1, Valor1, Coord2, Valor2)


	//Largo
	//(CoordObjetivo, Velocidad, Coord1, Valor1, Coord2, Valor2, Coord3, Valor3, Coord4, Valor4, EsquemaLargo)

	if(Rank != 0 && Rank != Procesos - 1){

		
		for(i = Ix; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
					wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

					wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
					wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

					uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

					wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
					wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

					wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
					wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior 
			j = 0;
			for(k = 1; k < NZ; k++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k-1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k-1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}

		}
	}	
	else if(Rank == 0){

		
		for(i = Ix + 1; i < Fx; i++){

			//Centro
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
					wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

					wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
					wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

					uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

					wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
					wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

					wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
					wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ; k++){
				
				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
			}


		}

		//Parte Izquierda
		i = 0;

		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){

				wW_pred = Wleft[WLEFT(i,j,k)];
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

				wW = Wleft[WLEFT(i,j,k)];
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ; k++){
				
				

				wW_pred = Wleft[WLEFT(i,j,k)];
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

				wW = Wleft[WLEFT(i,j,k)];
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

		}

		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){
				
				wW_pred = Wleft[WLEFT(,j,k)];
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = Wleft[WLEFT(i,j,k)];
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = 0.50*(Uleft[ULEFT(i,j,k)] + Uleft[ULEFT(i,j,k - 1)]);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
		}

	}
	else if (Rank == Procesos - 1){

		
		for(i = Ix; i < Fx - 1; i++){

			//Centro
			for(k = 1; k < NZ; k++){
				for(j = 1; j < NY - 1; j++){

					wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
					wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

					wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
					wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

					wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
					wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

					uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
					uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

					vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
					vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

					wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
					wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

					wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
					wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

					wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
					wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

					uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
					uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

					vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
					vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

					ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

				}
			}

			//Partes Inferior y Superior

			//Parte Inferior
			j = 0;

			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}


			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
			}


		}

		//Parte Derecha
		i = NX - 1;

		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Interpolacion(MESH.MU[GU(i+1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)]);

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[GU(i+1,j,k,0)], uE_pred, MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], MESH.MW[GW(i+2,j,k,0)], WFIELD[LW(i+2,j,k,0)], EsquemaLargo);

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wE_pred, MESH.MU[GU(i+1,j,k-2,2)], UFIELD[LU(i+1,j,k-2,0)], MESH.MU[GU(i+1,j,k-1,2)], UFIELD[LU(i+1,j,k-1,0)], MESH.MU[GU(i+1,j,k,2)], UFIELD[LU(i+1,j,k,0)], MESH.MU[GU(i+1,j,k+1,2)], UFIELD[LU(i+1,j,k+1,0)], EsquemaLargo);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

			}
		}

		//Partes Inferior y Superior

		//Parte Inferior
		j = 0;
		for(k = 1; k < NZ; k++){
				
				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Wright[WRIGHT(i,j,k)];

				wS_pred = Wbot[WBOT(i,j,k)];
				wN_pred = Interpolacion(MESH.MV[GV(i,j+1,k,1)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)]);

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS_pred = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = Wright[WRIGHT(i,j,k)];

				wS = Wbot[WBOT(i,j,k)];
				wN = ConvectiveScheme(MESH.MV[GV(i,j+1,k,1)], vN_pred, MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], MESH.MW[GW(i,j+2,k,1)], WFIELD[LW(i,j+2,k,0)], EsquemaLargo);

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS = 0.50*(Vbot[VBOT(i,j,k)] + Vbot[VBOT(i,j,k - 1)]);
				vN = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wN_pred, MESH.MV[GV(i,j+1,k-2,2)], VFIELD[LV(i,j+1,k-2,0)], MESH.MV[GV(i,j+1,k-1,2)], VFIELD[LV(i,j+1,k-1,0)], MESH.MV[GV(i,j+1,k,2)], VFIELD[LV(i,j+1,k,0)], MESH.MV[GV(i,j+1,k+1,0)], VFIELD[LV(i,j+1,k+1,0)], EsquemaLargo);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );
		}


		//Parte Superior
		j = NY - 1;
		for(k = 1; k < NZ; k++){

				wW_pred = Interpolacion(MESH.MU[GU(i,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)]);
				wE_pred = Wright[WRIGHT(i,j,k)];

				wS_pred = Interpolacion(MESH.MV[GV(i,j,k,1)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)]);
				wN_pred = Wtop[WTOP(i,j,k)];

				wH_pred = Interpolacion(MESH.MP[GP(i,j,k-1,2)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[GP(i,j,k,2)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)]);

				uW_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)]);
				uE_pred = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS_pred = Interpolacion(MESH.MW[GW(i,j,k,2)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)]);
				vN_pred = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				wW = ConvectiveScheme(MESH.MU[GU(i,j,k,0)], uW_pred, MESH.MW[GW(i-2,j,k,0)], WFIELD[LW(i-2,j,k,0)], MESH.MW[GW(i-1,j,k,0)], WFIELD[LW(i-1,j,k,0)], MESH.MW[GW(i,j,k,0)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i+1,j,k,0)], WFIELD[LW(i+1,j,k,0)], EsquemaLargo);
				wE = Wright[WRIGHT(i,j,k)];

				wS = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], vS_pred, MESH.MW[GW(i,j-2,k,1)], WFIELD[LW(i,j-2,k,0)], MESH.MW[GW(i,j-1,k,1)], WFIELD[LW(i,j-1,k,0)], MESH.MW[GW(i,j,k,1)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j+1,k,1)], WFIELD[LW(i,j+1,k,0)], EsquemaLargo);
				wN = Wtop[WTOP(i,j,k)];

				wH = ConvectiveScheme(MESH.MP[GP(i,j,k-1,2)], wH_pred, MESH.MW[GW(i,j,k-2,2)], WFIELD[LW(i,j,k-2,0)], MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[GP(i,j,k,2)], wT_pred, MESH.MW[GW(i,j,k-1,2)], WFIELD[LW(i,j,k-1,0)], MESH.MW[GW(i,j,k,2)], WFIELD[LW(i,j,k,0)], MESH.MW[GW(i,j,k+1,2)], WFIELD[LW(i,j,k+1,0)], MESH.MW[GW(i,j,k+2,2)], WFIELD[LW(i,j,k+2,0)], EsquemaLargo);

				uW = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wW_pred, MESH.MU[GU(i,j,k-2,2)], UFIELD[LU(i,j,k-2,0)], MESH.MU[GU(i,j,k-1,2)], UFIELD[LU(i,j,k-1,0)], MESH.MU[GU(i,j,k,2)], UFIELD[LU(i,j,k,0)], MESH.MU[GU(i,j,k+1,2)], UFIELD[LU(i,j,k+1,0)], EsquemaLargo);
				uE = 0.50*(Uright[URIGHT(i,j,k)] + Uright[URIGHT(i,j,k - 1)]);

				vS = ConvectiveScheme(MESH.MW[GW(i,j,k,2)], wS_pred, MESH.MV[GV(i,j,k-2,2)], VFIELD[LV(i,j,k-2,0)], MESH.MV[GV(i,j,k-1,2)], VFIELD[LV(i,j,k-1,0)], MESH.MV[GV(i,j,k,2)], VFIELD[LV(i,j,k,0)], MESH.MV[GV(i,j,k+1,0)], VFIELD[LV(i,j,k+1,0)], EsquemaLargo);
				vN = 0.50*(Vtop[VTOP(i,j,k)] + Vtop[VTOP(i,j,k - 1)]);

				ConvectiveW[LW(i,j,k,0)] = (1/MESH.VolMW[GW(i,j,k,0)])*(
											 - MESH.SupMW[GW(i,j,k,0)]*wW*uW
											 + MESH.SupMW[GW(i,j,k,1)]*wE*uE
											 - MESH.SupMW[GW(i,j,k,2)]*wS*vS
											 + MESH.SupMW[GW(i,j,k,3)]*wN*vN
											 - MESH.SupMW[GW(i,j,k,4)]*wH*wH
											 + MESH.SupMW[GW(i,j,k,5)]*wT*wT
											 );

		}

	}

}

//Cálculo de término de Boussinesq de la ecuación de Cantidad de Movimiento
void Solver::Get_BoussinesqTerm(Mesher MESH){
int i, j, k;
double Tv;

	//Velocidades V
	for(i = Ix; i < Fx; i++){
		for(k = 0; k < NZ; k++){
			for(j = 1; j < NY; j++){

				Tv = ConvectiveScheme(MESH.MV[GV(i,j,k,1)], VLPRES[LV(i,j,k,0)], MESH.MP[GP(i,j-2,k,1)], TLPRES[LP(i,j-2,k,0)], MESH.MP[GP(i,j-1,k,1)], TLPRES[LP(i,j-1,k,0)], MESH.MP[GP(i,j,k,1)], TLPRES[LP(i,j,k,0)], MESH.MP[GP(i,j+1,k,1)], TLPRES[LP(i,j+1,k,0)], EsquemaLargo);
					
				BoussinesqV[LV(i,j,k,0)] = gy*(1.0 - Beta*(Tv - To));

			}
		}
	}
	
}





//Resolución del sistema con Gauss-Seidel
void Solver::Get_GaussSeidel(ParPro MPI1){
int i, j, k;
MaxDiffGS = 2.0*ConvergenciaGS;

	while(MaxDiffGS >= ConvergenciaGS){

		MPI1.CommunicateDataLP(PLFUT, PLFUT, Ix, Fx);

		if(Rank != 0 && Rank != Procesos - 1){

			//Parte Central
			for(i = Ix; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}
			}
			
			//Parte Inferior
			j = 0;
			for(i = Ix; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte Superior
			j = NY - 1;
			for(i = Ix; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}
			
			//Parte Here
			k = 0;
			for(i = Ix; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
						PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte There
			k = NZ - 1;
			for(i = Ix; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Esquinas
			for(i = Ix; i < Fx; i++){

				//Esquina Abajo Here
				j = 0;
				k = 0;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
			}
			
		}
        else if(Rank == 0){
				
        	//Parte Central
			for(i = Ix + 1; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}
			}
			
			//Parte Inferior
			j = 0;
			for(i = Ix + 1; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte Superior
			j = NY - 1;
			for(i = Ix + 1; i < Fx; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte Here
			k = 0;
			for(i = Ix + 1; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
						PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte There
			k = NZ - 1;
			for(i = Ix + 1; i < Fx; i++){
				for(j = 1; j < NY - 1; j++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Esquinas
			for(i = Ix + 1; i < Fx; i++){

				//Esquina Abajo Here
				j = 0;
				k = 0;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
			}
			
			
			//Parte Izquierda
			i = 0;

			//Centro
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Esquina Abajo Here
			j = 0;
			k = 0;
			PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;
			PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
			//Esquina Abajo There
			j = 0;
			k = NZ - 1;
			PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;
			PLFUT[LP(i,j,k,0)] = (ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
		}
		else if(Rank == Procesos - 1){

			//Parte Central
			for(i = Ix; i < Fx - 1; i++){
				for(j = 1; j < NY - 1; j++){
					for(k = 1; k < NZ - 1; k++){
						PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
					}
				}
			}
			
			//Parte Inferior
			j = 0;
			for(i = Ix; i < Fx - 1; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte Superior
			j = NY - 1;
			for(i = Ix; i < Fx - 1; i++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte Here
			k = 0;
			for(i = Ix; i < Fx - 1; i++){
				for(j = 1; j < NY - 1; j++){
						PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte There
			k = NZ - 1;
			for(i = Ix; i < Fx - 1; i++){
				for(j = 1; j < NY - 1; j++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Esquinas
			for(i = Ix; i < Fx - 1; i++){

				//Esquina Abajo Here
				j = 0;
				k = 0;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Arriba Here
				j = NY - 1;
				k = 0;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Abajo There
				j = 0;
				k = NZ - 1;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
				//Esquina Arriba There
				j = NY - 1;
				k = NZ - 1;
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + ae[LA(i,j,k,0)]*PLFUT[LP(i+1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				
			}

			//Parte Derecha
			i = NX - 1;

			//Centro
			for(j = 1; j < NY - 1; j++){
				for(k = 1; k < NZ - 1; k++){
					PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
				}
			}

			//Parte Inferior
			j = 0;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Parte Superior
			j = NY - 1;
			for(k = 1; k < NZ - 1; k++){
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Parte Here
			k = 0;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Parte There
			k = NZ - 1;
			for(j = 1; j < NY - 1; j++){
				PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			}

			//Esquina Abajo Here
			j = 0;
			k = 0;
			PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
		
			//Esquina Arriba Here
			j = NY - 1;
			k = 0;
			PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,NZ - 1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,k+1,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
		
			//Esquina Abajo There
			j = 0;
			k = NZ - 1;

			PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + an[LA(i,j,k,0)]*PLFUT[LP(i,j+1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			
			//Esquina Arriba There
			j = NY - 1;
			k = NZ - 1;
			PLFUT[LP(i,j,k,0)] = (aw[LA(i,j,k,0)]*PLFUT[LP(i-1,j,k,0)] + as[LA(i,j,k,0)]*PLFUT[LP(i,j-1,k,0)] + ah[LA(i,j,k,0)]*PLFUT[LP(i,j,k-1,0)] + at[LA(i,j,k,0)]*PLFUT[LP(i,j,0,0)] + bp[LA(i,j,k,0)])/ap[LA(i,j,k,0)];
			
		}

		MaxDiffGS = 0.0;

		for(i = Ix; i < Fx; i++){
			for(k = 0; k < NZ; k++){
				for(j = 0; j < NY; j++){
					if(abs(PLFUT[LP(i,j,k,0)] - PLSUP[LP(i,j,k,0)]) >= MaxDiffGS){
						MaxDiffGS = abs(PLFUT[LP(i,j,k,0)] - PLSUP[LP(i,j,k,0)]);
						
					}
				}
			}
		}

		MPI_Allreduce(&MaxDiffGS, &MaxDiffGS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		//cout<<"MaxDiff: "<<MaxDiffGS<<endl;
		for(i = Ix; i < Fx; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					PLSUP[LP(i,j,k,0)] = PLFUT[LP(i,j,k,0)];	
				}
			}
		}

	}

}









			