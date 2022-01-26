#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <stdio.h>
#include <mpi.h>

using namespace std;

class Solver{
	private:

	public:
		Solver(Memory, ReadData, ParPro, Mesher, PostProcessing, string);

		

		

		//Nodos del dominio (Problema Square Cylinder)
		int NX1;
		int NX2;
		int NX3;

		int NY1;
		int NY2;
		int NY3;

		
		int StepsPantalla;
		int StepsFile;

		double ConvergenciaGS;
		double MaxDiffGS;

		double ConvergenciaGlobal;
		double MaxDiffGlobal;

		

		string DIRECTORIO;

		

		double Xcentroide;
		double Ycentroide;

		double Xcuadrado;
		double Ycuadrado;

		
		
		

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

		

		//Variables del Runge Kutta (Wikipedia)
		double c1, c2, c3, c4;
		double b1, b2, b3, b4;
		double a_21, a_22, a_23, a_24;
		double a_31, a_32, a_33, a_34;
		double a_41, a_42, a_43, a_44;

		//Matrices del Solver

		//Matrices globales de propiedades
		double *UGlobal; //Velocidad U Global
		double *VGlobal; //Velocidad V Global
		double *WGlobal; //Velocidad W Global
		double *PGlobal; //Presion P Global

		//Matrices locales de propiedades en los nodos de Presión

		//Presión
		double *PLFUT;
		double *PLSUP;

		//Velocidad U
		double *ULPRES;
		double *ULFUT;

		//Velocidad V
		double *VLPRES;
		double *VLFUT;

		//Velocidad W
		double *WLPRES;
		double *WLFUT;

		//Temperatura
		double *TLPAS;
		double *TLPRES;
		double *TLFUT;

		//Velocidad Predictoras
		double *PredU;
		double *PredV;
		double *PredW;

		double *K1_U;
		double *K1_V;
		double *K1_W;
		double *K1_T;

		double *K2_U;
		double *K2_V;
		double *K2_W;
		double *K2_T;

		double *K3_U;
		double *K3_V;
		double *K3_W;
		double *K3_T;

		double *K4_U;
		double *K4_V;
		double *K4_W;
		double *K4_T;

		double *UL_NEW;
		double *VL_NEW;
		double *WL_NEW;
		double *TL_NEW;

		double *ContribucionUpres;
		double *ContribucionUpast;

		double *ContribucionVpres;
		double *ContribucionVpast;

		double *ContribucionWpres;
		double *ContribucionWpast;

		//Matrices de los coeficientes de discretización de presión
		double *aw;
		double *ae;
		// Hacer una estructura para los coeficientes
		double *as;
		double *an;

		double *ah;
		double *at;

		double *ap;
		double *bp;

		//double *bpGlobal;
		
		//Matrices de las contribuciones de cada una de las ecuaciones
		double *ConvectiveU;
		double *ConvectiveV;
		double *ConvectiveW;

		double *DiffusiveU;
		double *DiffusiveV;
		double *DiffusiveW;

		double *ConvectiveT;
		double *DiffusiveT;

		double *TcontributionPast;
		double *TcontributionPres;

		double *BoussinesqU;
		double *BoussinesqV;
		double *BoussinesqW;

		//Matrices de condiciones de contorno
		
		//Velocidad U
		double *Uleft;
		double *Uright;
		double *Ubot;
		double *Utop;
		double *Uhere;
		double *Uthere;

		//Velocidad V
		double *Vleft;
		double *Vright;
		double *Vbot;
		double *Vtop;
		double *Vhere;
		double *Vthere;

		//Velocidad W
		double *Wleft;
		double *Wright;
		double *Wbot;
		double *Wtop;
		double *Where;
		double *Wthere;

		//Temperatura
		double *TLEFT;
		double *TRIGHT;
		double *TBOT;
		double *TTOP;
		double *There;
		double *Tthere;

		void AllocateMatrix(Memory);
		void Get_BoundaryConditions(Mesher);
		void Get_HaloVelocities();
		void Get_InitialConditions();
		void Get_StepTime(Mesher, ParPro);
		
		inline double ConvectiveScheme(double, double, double, double, double, double, double, double, double, double, string);
		
		void Get_PressureCoefficients(Mesher);

		void Get_DiffusiveU(Mesher, double*);
		void Get_DiffusiveV(Mesher, double*);
		void Get_DiffusiveW(Mesher, double*);
		void Get_DiffusiveT(Mesher, double*);

		void Get_ConvectiveU(Mesher, double*, double*, double*);
		void Get_ConvectiveV(Mesher, double*, double*, double*);
		void Get_ConvectiveW(Mesher, double*, double*, double*);
		void Get_ConvectiveT(Mesher, double*, double*, double*, double*);

		void Get_BoussinesqTerm(Mesher);

		void Get_Contributions(Mesher, ParPro);

		//void Get_RungeKuttaVelocities(Mesher, ParPro);
		//void Get_RungeKuttaTemperature(Mesher, ParPro);

		void Get_PredictorsDivergence(Mesher, ParPro);

		void Get_MaxDifGS(ParPro);
		void Get_GaussSeidel(ParPro);

		void Get_Velocities(Mesher, ParPro);

		void Get_Stop(ParPro MPI1);
		void Get_Update();

		void ExecuteSolver(Memory, ParPro, Mesher, PostProcessing);
		
};
