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

        //Matrices necesarias
        double* Test_Mesh;
        double* Test_MeshGlobal;

        
        // Mass Conservation
        double *Density;
        double *Density_Convective;

        // Structure for the properties values a the walls
        struct Walls
        {
            double *Wall_U;
            double *Wall_V;
            double *Wall_W;
        };

        // Structure for the boundary conditions
        struct Boundaries
        {
            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;
        };

        struct Walls Rho_W;
        struct Boundaries Rho_Boun;

        struct Walls U_W;
        struct Boundaries U_Boun;

        struct Walls V_W;
        struct Boundaries V_Boun;

        struct Walls W_W;
        struct Boundaries W_Boun;

        

		//Constructor de la clase
		CFD_Solver(Memory, ReadData, Parallel);
		
		//Metodos de la clase
        void AllocateMemory(Memory);
        inline double CS(double, double, double, double, double, double, double, double, double, double);      
		void Get_BoundaryConditions();
        void CommunicateVelocities(Parallel, double*, double*, double*);
        void ApplyBoundaries(Walls&, Boundaries&);
        void Get_WallsProperty_Scalar(Parallel, double*, double*, double*, double*, double*, Walls&, Boundaries&);
        void Get_WallsVelocities(double*, double*, double*, double*, Walls&, Boundaries&, Walls&, Boundaries&, Walls&, Boundaries&);
        void Get_ConvectiveTermScalar(Mesher, Walls&, double*);

};