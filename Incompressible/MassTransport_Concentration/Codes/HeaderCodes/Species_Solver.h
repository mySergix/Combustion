//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR SPECIES SOLVER CLASS                               //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>



using namespace std;

// Forward declaration of classes
class CFD_Solver;

class Species_Solver{
    private:

    public:

        Species_Solver(Memory, ReadData, Parallel);

        // Geometry Data
		double Xdominio;
		double Ydominio;
		double Zdominio;

        // Problem Data
		int NX;
		int NY; 
		int NZ;

        string EsquemaLargo;
		string EsquemaCorto;

        // Parameters for parallel computing
		int Rango;
		int Procesos;
        int Halo;
        int HP;

		int *Ix;
		int *Fx;
		
        // Physical Data
		double Rho;
		double Uref;
		double Reynolds;

		double Rayleigh;
		double Cp;
		double Prandtl;

		double gx;
		double gy;
		double gz;

		double DeltaT;
        double mu;

      

        

        // Class functions

            
};