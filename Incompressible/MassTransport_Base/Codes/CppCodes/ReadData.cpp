//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR READ INPUT DATA CLASS                                 //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"

using namespace std;

//Constructor del lector de datos
ReadData::ReadData(Memory M1){

	GeometryData = M1.AllocateDouble(7, 1, 1, 1); //Datos de la geometría del problema
	ProblemNumericalData = M1.AllocateInt(10, 1, 1, 1); //Datos numéricos del problema
	ProblemData = M1.AllocateDouble(5, 1, 1, 1); //Datos del problema
	ProblemPhysicalData = M1.AllocateDouble(11, 1, 1, 1); //Datos físicos sobre las condiciones del problema
	
}

void ReadData::ReadArrays(string FileName, int TotalData, double *Array){
int i = 0;
stringstream InitialName;
string FinalName;

	InitialName<<"../InputData/"<<FileName;
	FinalName=InitialName.str();

	ifstream Data(FinalName.c_str());

		if (Data){
        		string line;
        		while (getline(Data, line)){
        	 		istringstream iss(line);
					if(i < TotalData){
						if (iss >> Array[i]){ i++; }	
					}
        		 }
   	 	}
		
    	Data.close();

}

void ReadData::ReadInputs(){

	//Lectura datos en Arrays
	ReadArrays("GeometryData.txt", 7, GeometryData); //Input Datos Geometría del problema
	ReadArrays("ProblemPhysicalData.txt", 11, ProblemPhysicalData);	  //Input Datos físicos de las condiciones del problema

int i = 0;
string FileName;
stringstream InitialName;
string FinalName;

	FileName = "ProblemData.txt";

	InitialName<<"../InputData/"<<FileName;
	FinalName=InitialName.str();
	
	ifstream DatosProblema(FinalName.c_str());
	
		if (DatosProblema){
        		string line;
        		while (getline(DatosProblema, line)){
        	 		istringstream iss(line);
					if(i < 10){
						if (iss >> ProblemNumericalData[i]){ i++; }	
					}
					else if(i >= 10 && i < 15){
						if (iss >> ProblemData[i-10]){ i++; }	
					}
					else if(i == 15){
						if (iss >> ConvectiveScheme1){ i++; }	
					}
					else{
						if (iss >> ConvectiveScheme2){ i++; }	
					}
        		 }
   	 	}
			
    	DatosProblema.close();	

}

