//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR READ INPUT DATA CLASS                                 //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"

using namespace std;

//Constructor del lector de datos
ReadData::ReadData(Memory M1){

    GeometryData = M1.AllocateDouble(3, 1, 1, 1); // Geometry Data of the problem
	ProblemPhysicalData = M1.AllocateDouble(5, 1, 1, 1); // Physical Data of the problem
	MeshingOptionsData = M1.AllocateDouble(3, 1, 1, 1); //Datos del problema

	ProblemNumericalData = M1.AllocateInt(8, 1, 1, 1); //Datos numéricos del problema
   
}

void ReadData::ReadChemicalSpecies(Memory M1){
int i = 0;
int NumberSpecies = 0;

    // Get the total number of species present in the txt file
    std::string line;
    std::ifstream SpeciesFile("../InputData/ChemicalSpecies.txt");

    while (std::getline(SpeciesFile, line)){
        ++NumberSpecies;
    }

    NumberSpecies -= 3; //Initial Lines of the file

    ChemicalSpecies = M1.AllocateStringArray(NumberSpecies);

    SpeciesFile.clear();
    SpeciesFile.seekg(0,std::ios::beg);

    while (std::getline(SpeciesFile, line)){

        if (i >= 3){
            ChemicalSpecies[i - 3] = line;
        }
        i++;

    }

}

void ReadData::ReadDataFileToDoubleArray(string FileName, double* Array, int TotalData){
int i = 0;

stringstream InitialName;
string FinalName;

	InitialName<<"../InputData/"<<FileName;
	FinalName = InitialName.str();

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

void ReadData::ReadDataFileToIntArray(string FileName, int* Array, int TotalData){
int i = 0;
stringstream InitialName;
string FinalName;

	InitialName<<"../InputData/"<<FileName;
	FinalName = InitialName.str();

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

void ReadData::ReadInputs(Memory M1){

    ReadChemicalSpecies(M1); // Read Data from the species present file

    ReadDataFileToDoubleArray("GeometryData.txt", GeometryData, 3);
    ReadDataFileToDoubleArray("ProblemPhysicalData.txt", ProblemPhysicalData, 5);
    ReadDataFileToDoubleArray("MeshingOptionsData.txt", MeshingOptionsData, 3);

    ReadDataFileToIntArray("ProblemNumericalData.txt", ProblemNumericalData, 8);

}