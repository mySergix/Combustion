//------------------------------------------------------------------------------------------------//
//                     CPP FILE FOR SPECIES SOLVER CLASS DATA INPUT                               //
//------------------------------------------------------------------------------------------------//

// Function to read the name of the files of species data
void Species_Solver::Read_SpeciesName(string FileName){
int i = 0;
int lines = 0;
string linea;
stringstream InitialName;
string FinalName;

	InitialName<<"../InputData/SpeciesData/"<<FileName;
	FinalName = InitialName.str();

	ifstream Data(FinalName.c_str());

		if(Data){
			while(getline(Data, linea)){
				lines++;
			}
		}

		Data.close();

		SpeciesString = new string[lines - 2];

		ifstream Data2(FinalName.c_str());

		if (Data2){
			string line;
			while (getline(Data2, line)){
				istringstream iss(line);
				if(i >= 2 && i < lines){
					iss >> SpeciesString[i - 2];
				}
				i++;
			}
		}

		Data2.close();

}

// Function to read the data of each species
void Species_Solver::Read_SpeciesInformation(Species_Struct &Specie, string FileName){
int i = 0;
string linea;
stringstream InitialName;
string FinalName;

	InitialName<<"../InputData/SpeciesData/"<<FileName;
	FinalName = InitialName.str();

	ifstream Data(FinalName.c_str());

		if (Data){
        	string line;
			for (i = 0; i < 4; i++){
				getline(Data, line);
			}
			getline(Data, Specie.Name);
		}
		i = 0;
		if (Data){
        	string line;
        	while (getline(Data, line)){
        	 	istringstream iss(line);
				if(i == 0){
					if (iss >> Specie.Wmolar){ i++; }	
				}
				else if(i == 1){
					if (iss >> Specie.Epsilon){ i++; }	
				}
				else if(i == 2){
					if (iss >> Specie.sigma){ i++; }	
				}
				else if(i >= 3 && i < 11){
					if (iss >> Specie.mu_coeff[i-3]){ i++; }	
				}
				else if(i >= 11 && i < 19){
					if (iss >> Specie.lambda_coeff[i-11]){ i++; }	
				}
				else if(i >= 19 && i < 29){
					if (iss >> Specie.Cp_coeff[i-19]){ i++; }	
				}
				else if(i >= 29 && i < 41){
					if (iss >> Specie.h_coeff[i-29]){ i++; }	
				}
        	}
   	 	}
		
    Data.close();

}

// Function to read ALL of species data
void Species_Solver::Read_AllSpeciesData(){
int i;

    for (i = 0; i < N_Species; i++){
        Read_SpeciesInformation(Species[i], SpeciesString[i]);
    }

}