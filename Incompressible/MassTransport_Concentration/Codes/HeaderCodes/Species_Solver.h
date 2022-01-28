//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR SPECIES SOLVER CLASS                               //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#define N_Species 2

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

        double To;
		double Beta;
		double Difference;
		double Producto;

        double nu;

	    double D_AB;

		double K;
        double Co;
        double hfg;
        double T_water;

        double Schmidt;
        double MW_H2O;

        struct Species_Struct
        {
            double *C_Pres;
            double *C_Fut;

            double *ContributionPast;
            double *ContributionPres;

            double *Convective;
            double *Diffusive;

            double *D_ab;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            double *Global;
        };

        Species_Struct Species[N_Species];

        // Class functions

            // Memory Allocation
            void Allocate_StructSpecies(Memory, int);
            
            // Boundary Conditions
            void Get_InitialConditionsSpecies();
            void Get_InitialBoundaryConditionsSpecies();
            void Get_UpdateBoundaryConditionsSpecies(Mesher, double*);
            void Get_InitialHalosSpecies();
            void Get_UpdateHalosSpecies();

            // Diffusion
            void Get_DiffusionCoefficients(int);
            void Get_DiffusionSpecies(Mesher, int);
            void Get_MassConservationSpecies(int);

            // Convection
            inline double ConvectiveSchemeSpecies(double, double, double, double, double, double, double, double, double, double, string);
            void Get_ConvectionSpecies(Mesher, CFD_Solver, int);

            void Get_Concentration();
};