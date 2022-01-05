//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR SPECIES SOLVER CLASS                               //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

class Species_Solver{	

	public:

        // Variables de la clase

        // Geometry data
        double Xdominio;
        double Ydominio;
        double Zdominio;

        // Meshing data
        int NX;
        int NY;
        int NZ; 

        // Parallel computing data
        int Rango;
        int Procesos;
        int* Ix;
        int* Fx;

        int Halo;

        double *kB;

        //Matrices necesarias
        double* Test_Mesh;
        double* Test_MeshGlobal;

         // Structure for species data
        struct Species_Struct
        {
            string Name; // Name of the species
            double *Wmolar; // Molar weight of the species
            double *Epsilon; // Characteristic Lennard-Jones energy

            // Mass fractions
            double *Y_Past;
            double *Y_Pres;
            double *Y_Fut;

            double *Y_Wall_U;
            double *Y_Wall_V;
            double *Y_Wall_W;

            // Molar fraction
            double *X;

            // Diffusion velocities
            double *U_Diff;
            double *V_Diff;
            double *W_Diff;

            // Equation terms
            double *ConvectiveTerm;

            double *ContributionPres;
            double *ContributionPast;
        };

        struct Species_Struct Species[N_Species];

		//Constructor de la clase
		Species_Solver(Memory M1, ReadData R1, Parallel P1);
		
		//Metodos de la clase

        // Binary Diffusion Coefficient Models
        double Get_BinaryDiff_ChampanEnskog(int, double, double, int, int, int);
        double Get_BinaryDiff_WilkeLee(int, double, double, int, int, int);

        // Diffusion Models
        void Get_DiffusionCoefficient_FickModel(int, CFD_Solver);

        void Get_WallsDiffusionVelocities(int, Mesher);

        void Get_SpeciesDiffusion(Mesher, CFD_Solver, int);
        void Get_SpeciesConvection(Mesher, CFD_Solver, int);
        void Get_StepContribution_Species(CFD_Solver, int);
        void Get_TemporalIntegration_Species(int);
        void Get_Update(int);
        void Get_MolarFraction_X();

        // JANAF Calculations
        double JANAF_CpHeat(double, int, int, int);
        double JANAF_AbsEnthalpy_Specie(int, double);
        double JANAF_AbsEnthalpy_Specie_Mix(double, int, int, int);
        double JANAF_DynViscosity(double, int, int, int);
        double JANAF_ThermalCond(double, int, int, int);

};