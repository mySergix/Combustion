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

		string ConvectiveScheme1;
		string ConvectiveScheme2;
		
		double *GeometryData; //Datos de la geometría del problema
		double *ProblemPhysicalData; //Datos físicos sobre las condiciones del problema
		double *ProblemData; //Datos del problema
		int *ProblemNumericalData; //Datos numéricos del problema
		
		//Constructor de la clase
		ReadData(Memory);

		void ReadInputs(); //Lector datos en ficheros
		void ReadArrays(string, int , double*);

};
