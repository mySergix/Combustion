//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MESH CREATION CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"

using namespace std;

//Constructor del mallador
Mesher::Mesher(Memory M1, ReadData R1, Parallel P1){
		
	//Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	NX = R1.ProblemNumericalData[2];
	NY = R1.ProblemNumericalData[3];
	NZ = R1.ProblemNumericalData[4];
	Halo = 2;
	HP = 2;

	OptionX = R1.ProblemNumericalData[5];
	OptionY = R1.ProblemNumericalData[6];
	OptionZ = R1.ProblemNumericalData[7];

	SFX = R1.ProblemData[0];
	SFY = R1.ProblemData[1];
	SFZ = R1.ProblemData[2];

	//Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	//Datos necesarios para computación paralela
	Rango = P1.Rango;
	Procesos = P1.Procesos;
	Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

    for (int i = 0; i < Procesos; i++){
        Ix[i] = P1.Ix[i];
        Fx[i] = P1.Fx[i];
    }

}

#include "Matrix_Index.cpp"
#include "Mesher_Memory.cpp"
#include "Mesher_NodalCoordinates.cpp"

//Cálculo de las distancias entre nodos en cada una de las matrices
void Mesher::Get_Deltas(){
int i, j, k;

	// Collocated mesh distances
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				DeltasMP[LP(i,j,k,0)] = MU[LU(i+1,j,k,0)] - MU[LU(i,j,k,0)]; // Delta X
				DeltasMP[LP(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
				DeltasMP[LP(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
			}
		}
	}

	// Staggered U mesh distances
	if (Rango != 0 && Rango != Procesos - 1){
		for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					DeltasMU[LU(i,j,k,0)] = MP[LP(i,j,k,0)] - MP[LP(i-1,j,k,0)]; // Delta X
					DeltasMU[LU(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
					DeltasMU[LU(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
				}
			}
		}
	}
	else if (Rango == 0){ // Left
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){

				DeltasMU[LU(0,j,k,0)] = MP[LP(0,j,k,0)] - MU[LU(0,j,k,0)]; // Delta X Parte Izquierda
				DeltasMU[LU(0,j,k,1)] = MV[LV(0,j+1,k,1)] - MV[LV(0,j,k,1)]; // Delta Y Parte Izquierda
				DeltasMU[LU(0,j,k,2)] = MW[LW(0,j,k+1,2)] - MW[LW(0,j,k,2)]; // Delta Z Parte Izquierda

				for(i = Ix[Rango] + 1; i < Fx[Rango] + 1; i++){
					DeltasMU[LU(i,j,k,0)] = MP[LP(i,j,k,0)] - MP[LP(i-1,j,k,0)]; // Delta X
					DeltasMU[LU(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
					DeltasMU[LU(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
				}
			}
		}
	}
	else if(Rango == Procesos - 1){ // Right
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){

				DeltasMU[LU(NX,j,k,0)] = MU[LU(NX,j,k,0)] - MP[LP(NX-1,j,k,0)]; //Delta X Parte Derecha
				DeltasMU[LU(NX,j,k,1)] = MV[LV(NX-1,j+1,k,1)] - MV[LV(NX-1,j,k,1)]; //Delta Y Parte Derecha
				DeltasMU[LU(NX,j,k,2)] = MW[LW(NX-1,j,k+1,2)] - MW[LW(NX-1,j,k,2)]; //Delta Z Parte Derecha

				for(i = Ix[Rango]; i < Fx[Rango]; i++){
					DeltasMU[LU(i,j,k,0)] = MP[LP(i,j,k,0)] - MP[LP(i-1,j,k,0)]; // Delta X
					DeltasMU[LU(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
					DeltasMU[LU(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
				}
			}
		}
	}

	// Staggered V mesh distances
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(k = 0; k < NZ; k++){

			DeltasMV[LV(i,0,k,0)] = MU[LU(i+1,0,k,0)] - MU[LU(i,0,k,0)]; // Delta X Parte Abajo
			DeltasMV[LV(i,NY,k,0)] = MU[LU(i+1,NY-1,k,0)] - MU[LU(i,NY-1,k,0)]; // Delta X Parte Arriba

			DeltasMV[LV(i,0,k,1)] = MP[LP(i,0,k,1)] - MV[LV(i,0,k,1)]; // Delta Y Parte Abajo
			DeltasMV[LV(i,NY,k,1)] = MV[LV(i,NY,k,1)] - MP[LP(i,NY-1,k,1)]; // Delta Y Parte Arriba

			DeltasMV[LV(i,0,k,2)] = MW[LW(i,0,k+1,2)] - MW[LW(i,0,k,2)]; // Delta Z Parte Abajo
			DeltasMV[LV(i,NY,k,2)] = MW[LW(i,NY-1,k+1,2)] - MW[LW(i,NY-1,k,2)]; // Delta Z Parte Arriba

			for(j = 1; j < NY; j++){
				DeltasMV[LV(i,j,k,0)] = MU[LU(i+1,j,k,0)] - MU[LU(i,j,k,0)]; //Delta X
				DeltasMV[LV(i,j,k,1)] = MP[LP(i,j,k,1)] - MP[LP(i,j-1,k,1)]; //Delta Y
				DeltasMV[LV(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; //Delta Z
			}

		}
	}

	// Staggered W mesh distances
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){

			DeltasMW[LW(i,j,0,0)] = MU[LU(i+1,j,0,0)] - MU[LU(i,j,0,0)]; // Delta X Parte Here
			DeltasMW[LW(i,j,NZ,0)] = MU[LU(i+1,j,NZ-1,0)] - MU[LU(i,j,NZ-1,0)]; // Delta X Parte There

			DeltasMW[LW(i,j,0,1)] = MV[LV(i,j+1,0,1)] - MV[LV(i,j,0,1)]; // Delta Y Parte Here
			DeltasMW[LW(i,j,NZ,1)] = MV[LV(i,j+1,NZ-1,1)] - MV[LV(i,j,NZ-1,1)]; // Delta Y Parte There

			DeltasMW[LW(i,j,0,2)] = MP[LP(i,j,0,2)] - MW[LW(i,j,0,2)]; // Delta Z Parte Here
			DeltasMW[LW(i,j,NZ,2)] = MW[LW(i,j,NZ,2)] - MP[LP(i,j,NZ-1,2)]; // Delta Z Parte There

			for(k = 1; k < NZ; k++){
				DeltasMW[LW(i,j,k,0)] = MU[LU(i+1,j,k,0)] - MU[LU(i,j,k,0)]; //Delta X
				DeltasMW[LW(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; //Delta Y
				DeltasMW[LW(i,j,k,2)] = MP[LP(i,j,k,2)] - MP[LP(i,j,k-1,2)]; //Delta Z
			}
		}
	}

}

//Cálculo de las superficies de cada uno de los volúmenes de control
void Mesher::Get_Surfaces(){
int i, j, k;

	// Collocated mesh surfaces
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				SupMP[LP(i,j,k,0)] = DeltasMU[LU(i,j,k,1)]*DeltasMU[LU(i,j,k,2)]; // Surface X
				SupMP[LP(i,j,k,1)] = DeltasMV[LV(i,j,k,0)]*DeltasMV[LV(i,j,k,2)]; // Surface Y
				SupMP[LP(i,j,k,2)] = DeltasMW[LW(i,j,k,0)]*DeltasMW[LW(i,j,k,1)]; // Surface Z
			}
		}
	}

	// Staggered U mesh surfaces
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				SupMU[LU(i,j,k,0)] = DeltasMU[LU(i,j,k,1)]*DeltasMU[LU(i,j,k,2)]; // Surface X
				SupMU[LU(i,j,k,1)] = DeltasMU[LU(i,j,k,0)]*DeltasMU[LU(i,j,k,2)]; // Surface Y 
				SupMU[LU(i,j,k,2)] = DeltasMU[LU(i,j,k,0)]*DeltasMU[LU(i,j,k,1)]; // Surface Z
			}
		}
	}	

	// Staggered V mesh surfaces
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				SupMV[LV(i,j,k,0)] = DeltasMV[LV(i,j,k,1)]*DeltasMV[LV(i,j,k,2)]; // Surface X
				SupMV[LV(i,j,k,1)] = DeltasMV[LV(i,j,k,0)]*DeltasMV[LV(i,j,k,2)]; // Surface Y 
				SupMV[LV(i,j,k,2)] = DeltasMV[LV(i,j,k,0)]*DeltasMV[LV(i,j,k,1)]; // Surface Z
			}
		}
	}

	// Staggered W mesh surfaces
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){		
			for(k = 0; k < NZ + 1; k++){
				SupMW[LW(i,j,k,0)] = DeltasMW[LW(i,j,k,1)]*DeltasMW[LW(i,j,k,2)]; // Surface X
				SupMW[LW(i,j,k,1)] = DeltasMW[LW(i,j,k,0)]*DeltasMW[LW(i,j,k,2)]; // Surface Y 
				SupMW[LW(i,j,k,2)] = DeltasMW[LW(i,j,k,0)]*DeltasMW[LW(i,j,k,1)]; // Surface Z
			}
		}
	}

}

//Cálculo de los volúmenes de control de cada volúmen
void Mesher::Get_Volumes(){
int i, j, k;

	// Collocated mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				VolMP[LP(i,j,k,0)] = DeltasMP[LP(i,j,k,0)]*SupMP[LP(i,j,k,0)];
			}
		}
	}

	// Staggered U mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				VolMU[LU(i,j,k,0)] = DeltasMU[LU(i,j,k,0)]*SupMU[LU(i,j,k,0)];
			}
		}
	}

	// Staggered V mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){
				VolMV[LV(i,j,k,0)] = DeltasMV[LV(i,j,k,0)]*SupMV[LV(i,j,k,0)];
			}
		}
	}

	// Staggered W mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ + 1; k++){
				VolMW[LW(i,j,k,0)] = DeltasMW[LW(i,j,k,2)]*SupMW[LW(i,j,k,2)];
			}
		}
	}

}

//Ejecutar todos los procesos del mallador
void Mesher::ExecuteMesher(Memory M1){
	
	Allocate_MesherMemory(M1); // Memory Allocation
	Get_Meshes(); // Creación de todas las mallas
	Get_Deltas(); // Cálculo de las distancias entre nodos en cada una de las matrices
	Get_Surfaces(); // Cálculo de las superficies de cada uno de los volúmenes de control
	Get_Volumes(); // Cálculo de los volúmenes de control de cada volúmen

	if(Rango == 0){	
	//	MallaVTK3D("ParaviewResults/MeshResults/", "MallaP", "MallaMP", MP, NX, NY, NZ);
		cout<<"Mesh created."<<endl;
	}

}

