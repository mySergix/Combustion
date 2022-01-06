//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR CFD SOLVER CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

class CFD_Solver{	

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

        double Beta;

        double DeltaT;
        double MaxDiff;
        double GlobalConvergence;

        double c;
        double R_ideal;

        // Structure for the density, velocities and enthalpy maps
        struct Property_Struct
        {
            double *Past;
            double *Pres;
            double *Fut;

            double *ContributionPast;
            double *ContributionPres;

            double *Convective;
            double *Diffusive;

            double *Wall_U;
            double *Wall_V;
            double *Wall_W;

            double *Wall_U_Diff_X;
            double *Wall_U_Diff_Y;
            double *Wall_U_Diff_Z;

            double *Wall_V_Diff_X;
            double *Wall_V_Diff_Y;
            double *Wall_V_Diff_Z;

            double *Wall_W_Diff_X;
            double *Wall_W_Diff_Y;
            double *Wall_W_Diff_Z;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            double Gravity;
        };

        // Structure for the pressure field
        struct Pressure_Struct
        {
            double *Pres;

            double *Gradient_X;
            double *Gradient_Y;
            double *Gradient_Z;
        };

        // Structure for the viscous stresses
        struct Stresses_Struct
        {
            double *mu_Visc;

            double *Divergence;

            double *Tau_xx_X;
            double *Tau_yy_Y;
            double *Tau_zz_Z;

            double *Tau_yx_Y;
            double *Tau_zx_Z;

            double *Tau_xy_X;
            double *Tau_zy_Z;

            double *Tau_yz_Y;
            double *Tau_xz_X;

        };

        // Structure for the global matrix of core 0
        struct Global
        {
            double *Density;
            double *U;
            double *V;
            double *W;
        };

        // Aditional Energy Equation terms
        double *Lambda;
        double *FourierDiffusion;
        double *EnthalpyDiffusion;

        // Structures Declaration
        struct Property_Struct Density;

        struct Property_Struct U;
        struct Property_Struct V;
        struct Property_Struct W;

        struct Presion_Struct Pressure;
        
        struct Global GlobalMatrix;

        struct Stresses_Struct Stress;

		//Constructor de la clase
		CFD_Solver(Memory, ReadData, Parallel);
		
		//Metodos de la clase

            // Memory Allocation
            void Allocate_Struct_MapFields(Memory, Property_Struct&);
            void Allocate_Struct_VelocityGradients(Memory, Property_Struct&);
            void Allocate_Struct_BoundaryConditions(Memory, Property_Struct&);
            void Allocate_Struct_Contributions(Memory, Property_Struct&);
            void Allocate_Struct_VelocityEqsTerms(Memory, Property_Struct&);
            void Allocate_Struct_Pressure(Memory);
            void Allocate_Struct_Stresses(Memory);
            void Allocate_Struct_EnergyEqTerms(Memory);
            void Allocate_Struct_Global(Memory);

        inline double CS(double, double, double, double, double, double, double, double, double, double);   
        void Set_InitialValues();
		void Get_BoundaryConditions();
        void Update_BoundaryConditions(Mesher);
        void ApplyBoundaries(Property_Struct&);
        void Get_PeriodicConditions(Mesher, Property_Struct&);

        void Get_TimeStep(Mesher);

        void CommunicateVelocities(Parallel, double*, double*, double*);
        
        void Get_WallsValue_Property(Parallel, Mesher, Property_Struct&);

        void Get_ConvectiveTerm_Density(Mesher);
        void Get_ConvectiveTerm_Property(Mesher, Property_Struct&);
        void Get_Divergence(Mesher);
        void Get_VelocityGradients(Mesher, Property_Struct&);

        void Get_DynamicViscosity(Species_Solver);
        void Get_ViscousStresses(Mesher);
        void Get_DiffusiveTerm_Velocity(Mesher, Property_Struct&, double*, double*, double*);

        void Get_Pressure();
        void Get_PressureGradient(Mesher);

        void Get_StepContribution_Density(Property_Struct&);
        void Get_StepContribution_Velocity(Property_Struct&, double*);
        
        void Get_TemporalIntegration_Property(Property_Struct&);

        void Get_MaximumDifference(Property_Struct&, double&);
        void Get_ConvergenceCriteria();

        void UpdateField(Property_Struct&);

        void RunSolver(Memory, Parallel, Mesher, PostProcess);

};