//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR MESH CREATION CLASS                                //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

using namespace std;

class Mesher{
	private:

	public:
		//Constructor de la clase
		Mesher(Memory, ReadData, Parallel);

		// Datos Numéricos del problema
		int Problema;

		int ProgresionLineal;

		int NX;
		int NY; 
		int NZ;

		int OptionX;
		int OptionY;
		int OptionZ;

		double SFX;
		double SFY;
		double SFZ;

		// Datos Geométricos del problma
		double Xdominio;
		double Ydominio;
		double Zdominio;

		// Parámetros de computación paralela
		int Rango;
		int Procesos;
		int* Ix;
        int* Fx;

		int HP;
		int Halo;

		// Meshes of the simulation
		double *MP; // Collocated mesh
		double *MU; // Staggered U mesh
		double *MV; // Staggered V mesh
		double *MW; // Staggered W mesh

		// Meshes nodal distances 
		double *DeltasMP; // Collocated mesh distances
		double *DeltasMU; // Staggered U mesh distances
		double *DeltasMV; // Staggered V mesh distances
		double *DeltasMW; // Staggered W mesh distances

		// Meshes surfaces of the CV
		double *SupMP; // Collocated mesh surfaces
		double *SupMU; // Staggered U mesh distances
		double *SupMV; // Staggered V mesh distances
		double *SupMW; // Staggered W mesh distances

		// Volumes of the CV of the meshes
		double *VolMP; // Collocated mesh CV Volumes
		double *VolMU; // Staggered U mesh CV Volumes
		double *VolMV; // Staggered V mesh CV Volumes
		double *VolMW; // Staggered W mesh CV Volumes

		void Allocate_MesherMemory(Memory); //Alojamiento de memoria para cada matriz
		void Delete_MesherMemory(); // Borrado de toda la memoria

		void Get_Meshes(); //Creación de todas las mallas
		void Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
		void Get_Surfaces(); //Cálculo de las superficies de cada uno de los volúmenes de control
		void Get_Volumes(); //Cálculo de los volúmenes de control de cada volúmen
		
		//void MallaVTK3D(string, string, string, double*, int, int, int); //Pasar todos los resultados escalares en 2D a un archivo VTK
		
		void ExecuteMesher(Memory); //Ejecutar todos los procesos del mallador

};
		
