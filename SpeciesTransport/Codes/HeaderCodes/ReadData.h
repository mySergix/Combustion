//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR READ INPUT DATA CLASS                              //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h>
#include <string>

using namespace std;

class ReadData{
	private:
		
		
	public:

        string *ChemicalSpecies; // Chemical Species present in the reaction

		double *GeometryData; // Geometry Data of the problem
		double *ProblemPhysicalData; // Physical Data of the problem
		double *MeshingOptionsData; //Datos del problema

		int *ProblemNumericalData; //Datos numéricos del problema
		
		//Constructor de la clase
		ReadData(Memory);

        void ReadChemicalSpecies(Memory);

        void ReadDataFileToDoubleArray(string, double*, int);
        void ReadDataFileToIntArray(string, int*, int);

		void ReadInputs(Memory); // Lector datos en ficheros
		//void ReadArrays(string, int , double*);

		
};
