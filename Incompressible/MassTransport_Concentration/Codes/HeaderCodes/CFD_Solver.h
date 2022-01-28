//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR CFD SOLVER CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

// Forward declaration of classes
class Species_Solver;

class CFD_Solver{
    private:

    public:

        CFD_Solver(Memory, ReadData, Parallel);

        // Problem type
		int Problema;

        // Geometry Data
		double Xdominio;
		double Ydominio;
		double Zdominio;

        double Dh;
        double gamma_geometry;

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
        double nu;
        double mu;

        double Schmidt;
        double D_AB;
        double w_av;

        double Co;
        double To;
		double Beta;
		double Beta_Concentration;

		double K;
        double T_water;

        double ConvergenciaGS;
		double MaxDiffGS;

		double ConvergenciaGlobal;
		double MaxDiffGlobal;

        struct Velocity_Struct
        {
            double *Predictor;

            double *Pres;
            double *Fut;

            double *ContributionPast;
            double *ContributionPres;

            double *Convective;
            double *Diffusive;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            double *Boussinesq;
            double *Boussinesq_C;

            double Gravity;
        };

        struct Energy_Struct
        {
            double *Pres;
            double *Fut;

            double *ContributionPast;
            double *ContributionPres;

            double *Convective;
            double *Diffusive;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;
        };  

        struct Poisson_Coeffs
        {
            double *aw;
            double *ae;

            double *as;
            double *an;

            double *ah;
            double *at;

            double *ap;
            double *bp;
        };

        struct Pressure_Struct
        {
            double *Pres;
            double *Sup;
        };

        struct Global_Struct
        {
            double *P;
            double *U;
            double *V;
            double *W;

            double *T;
        };

        struct Velocity_Struct U;
        struct Velocity_Struct V;
        struct Velocity_Struct W;

        struct Pressure_Struct P;

        struct Energy_Struct T;
        struct Poisson_Coeffs A;
        
        struct Global_Struct Global;

        // Class functions

            // Memory Allocation
            void Allocate_PoissonCoeffsMemory(Memory);
            void Allocate_VelocitiesPartMemory(Memory, Velocity_Struct&, int, int, int);
            void Allocate_VelocitiesBoundaryConditionsMemory(Memory, Velocity_Struct&, int, int, int);
            void Allocate_PressureMemory(Memory);
            void Allocate_VelocitiesMemory(Memory);
            void Allocate_EnergyMemory(Memory);
            void Allocate_GlobalMemory(Memory);
            void Delete_VelocityMemory(Velocity_Struct&);

            // Utilities
            void Get_InitialConditions();
            void Get_StepTime(Mesher, Parallel);
            inline double ConvectiveScheme(double, double, double, double, double, double, double, double, double, double, string);
            void Get_ContributionsPredictors();
            void Get_PredictorsDivergence(Mesher);
            void Get_Velocities(Mesher, Parallel);
            void Get_Stop();
            void Get_Update();
            void Get_ConcentrationBuoyancy(Mesher, Species_Solver, int);

            // Boundary Conditions
            void Get_InitialBoundaryConditions(Mesher);
            void Get_UpdateBoundaryConditions(Mesher);
            void Get_StaticHalos();
            void Get_UpdateHalos();

            // Poisson Coefficients
            void Get_PoissonCoefficients(Mesher);
            void Get_GaussSeidel(Parallel);
            void PrintTxt();

            // Momentum Diffusion
            void Get_DiffusionU(Mesher);
            void Get_DiffusionV(Mesher);
            void Get_DiffusionW(Mesher);

            // Momentum Convection
            void Get_ConvectionU(Mesher);
            void Get_ConvectionV(Mesher);
            void Get_ConvectionW(Mesher);

            // Energy Equation
            void Get_DiffusionEnergy(Mesher);
            void Get_ConvectionEnergy(Mesher);
            void Get_BoussinesqV(Mesher);
            void Get_Temperature();

            // Run Solver
            void RunSolver(Memory, Parallel, Mesher, PostProcessing, Species_Solver);
            
};