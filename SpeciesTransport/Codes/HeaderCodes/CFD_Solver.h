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

        // Structure for the fields time steps properties
        struct Property
        {
            double *Past;
            double *Pres;
            double *Fut;

            double *ContributionPast;
            double *ContributionPres;

            double *Convective;

            double *Wall_U;
            double *Wall_V;
            double *Wall_W;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            double Gravity;
        };

        // Structure for the pressure field
        struct Presion
        {
            double *Pres;

            double *Gradient_X;
            double *Gradient_Y;
            double *Gradient_Z;
        };

        // Structure for the global matrix of core 0
        struct Global
        {
            double *Density;
            double *U;
            double *V;
            double *W;
        };

        struct Property Density;
        struct Property U;
        struct Property V;
        struct Property W;

        struct Presion Pressure;
        
        struct Global GlobalMatrix;

		//Constructor de la clase
		CFD_Solver(Memory, ReadData, Parallel);
		
		//Metodos de la clase
        void AllocateMemory(Memory);
        void Allocate_StructureMemory(Memory, Property&);

        inline double CS(double, double, double, double, double, double, double, double, double, double);   
        void Set_InitialValues();
		void Get_BoundaryConditions();
        void Update_BoundaryConditions();
        void Get_TimeStep(Mesher);

        void CommunicateVelocities(Parallel, double*, double*, double*);
        void ApplyBoundaries(Property&);
        void Get_WallsValue_Scalar(Parallel, double*, Property&);
        void Get_WallsVelocities(Mesher);
        void Get_ConvectiveTerm(Mesher, Property&);

        void Get_Pressure();
        void Get_PressureGradient(Mesher);
        void Get_StepContribution_Density(Property&);
        void Get_StepContribution_Velocity(Property&, double*);
        
        void Get_TemporalIntegration(Property&);

        void Get_MaximumDifference(Property&, double&);
        void Get_ConvergenceCriteria();

        void UpdateField(Property&);

        void RunSolver(Memory, Parallel, Mesher, PostProcess);

};