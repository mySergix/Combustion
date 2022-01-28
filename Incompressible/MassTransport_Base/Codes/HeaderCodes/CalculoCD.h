#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <stdio.h>
#include <mpi.h>

using namespace std;

class CalculoCD{
	private:

	public:
		CalculoCD(ReadData, ParPro, Mesher, Solver);

		//Datos Numéricos del problema
		int Problema;

		//Nodos del dominio
		int NX;
		int NY; 
		int NZ;

		string EsquemaLargo;
		string EsquemaCorto;

		//Parámetros de computación paralela
		int Rank;
		int Procesos;
		int Ix;
		int Fx;
		int Halo;

		int HP;
		
		//Datos Físicos del Problema
		double Rho;
		double Uref;
		double Reynolds;

		double Rayleigh;
		double Cp;
		double Prandtl;

		double gx;
		double gy;
		double gz;

		double Tleft;
		double Tright;

		double Tbot;
		double Ttop;

		double To;
		double Beta;
		double Difference;
		double Producto;

		double K;
		double mu;
	
		inline double ConvectiveScheme(double, double, double, double, double, double, double, double, double, double, string);
		
		void Get_DiffusiveU(Mesher, Solver, double&, double*);
		void Get_DiffusiveV(Mesher, Solver, double&, double*);
		void Get_DiffusiveW(Mesher, Solver, double&, double*);
		void Get_DiffusiveT(Mesher, Solver, double&, double*);

		void Get_ConvectiveU(Mesher, double*, double*, double*);
		void Get_ConvectiveV(Mesher, double*, double*, double*);
		void Get_ConvectiveW(Mesher, double*, double*, double*);
		void Get_ConvectiveT(Mesher, double*, double*, double*, double*);

		void Get_BoussinesqTerm(Mesher);
		
};