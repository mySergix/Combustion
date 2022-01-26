#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

#define PI 3.141592653589793

class PostProcessing{	

	private:

	public:
		//Constructor de la clase
		PostProcessing(Memory, ReadData, Mesher, string);
		
		//Datos de la clase
		int Problema;

    	int NX;
    	int NY;
    	int NZ;

    	double Reynolds;
    	double Rayleigh;
    	double Prandtl;

   		double Tleft; 
		double Tright;  

    	//Datos Geométricos del problema
		double Xdominio;
		double Ydominio;
		double Zdominio;

		double Xcentroide;
		double Ycentroide;

		double Xcuadrado;
		double Ycuadrado;

		int Halo;
		int HP;

		int HaloPressure;
		int HaloU;
		int HaloV; 
		int HaloW;
		
		string DIRECTORIO;
		
		//Metodos de la clase

		void LocalEscalarVTK3D(string, string, string, double*, double*, int, int, int, int, int);
		void LocalEscalarCoefficientsVTK3D(string, string, string, double*, double*, int, int, int, int);
		void LocalVectorialVTK3D(Mesher, string, string, string, double*, double*, double*, double*, int, int, int, int, int);

		void EscalarVTK3D(string, string, string, double*, double*, int, int, int);
		void VectorialVTK3D(Mesher, string, string, string, double*, double*, double*, double*, int, int, int);

		void Get_DrivenResults(Mesher, double, int, int, double*, double*);
		void Get_DifferentiallyResults(Mesher, double, int, int, double*, double*, double*);

};
