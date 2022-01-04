//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR JANAF PROPERTIES CALCULATIONS                         //
//------------------------------------------------------------------------------------------------//

// Cuidado, solo esta programado para 1 rango de temperaturas segun las tablas JANAF (200 - 1000), (1000 - 5000) K
// Falta de programar la evaluaccion de temperatura y decidir qué rango usar

// Function to calculate the Heat Capacity coefficient of the mixture in a control volume
double Species_Solver::JANAF_CpHeat(double T, int i, int j, int k){
int i;
double Cp_specie;
double Cp_global = 0.0;

    for (i = 0; i < N_Species; i++){
        Cp_specie = R_ideal * (Species[SP].Cp_coeff[0] + Species[SP].Cp_coeff[1] * T + Species[SP].Cp_coeff[2] * pow(T, 2) + Species[SP].Cp_coeff[3] * pow(T, 3) + Species[SP].Cp_coeff[4] * pow(T, 4));
        Cp_global += Species[SP].Y_Pres[LM(i,j,k,0)] * Cp_specie;
    }

    return Cp_global;

}

// Function to calculate the enthalpy of a specie for a certain Temperature
double Species_Solver::JANAF_AbsEnthalpy_Specie(int SP, double T){

    return R_ideal * T * (Species[SP].h_coeff[0] + (Species[SP].h_coeff[1]/2) * T + (Species[SP].h_coeff[2]/3) * pow(T, 2) + (Species[SP].h_coeff[3]/4) * pow(T, 3) + (Species[SP].h_coeff[4]/5) * pow(T, 4) + (Species[SP].h_coeff[1]/T));

}

// Function to calculate the absolute enthalpy of the mixture in a control volume
double Species_Solver::JANAF_AbsEnthalpy_Specie_Mix(double T, int i, int j, int k){
int i;
double AbsEnth;
double h = 0.0;

    for (i = 0; i < N_Species; i++){
        AbsEnth = R_ideal * T * (Species[SP].h_coeff[0] + (Species[SP].h_coeff[1]/2) * T + (Species[SP].h_coeff[2]/3) * pow(T, 2) + (Species[SP].h_coeff[3]/4) * pow(T, 3) + (Species[SP].h_coeff[4]/5) * pow(T, 4) + (Species[SP].h_coeff[1]/T));
        h += Species[SP].Y_Pres[LM(i,j,k,0)] * AbsEnth;
    }

    return h;

}

// Function to calculate the dynamic viscosity of the mixture in a control volume
double Species_Solver::JANAF_DynViscosity(double T, int i, int j, int k){
int i;
double DynVisc;
double mu = 0.0;

    for (i = 0; i < N_Species; i++){
        DynVisc = exp(Species[SP].mu_coeff[0] + Species[SP].mu_coeff[1] * log(T) + Species[SP].mu_coeff[2] * pow(log(T), 2) + Species[SP].mu_coeff[3] * pow(log(T), 3));
        mu += Species[SP].Y_Pres[LM(i,j,k,0)] * DynVisc;
    }

    return mu;

}

// Function to calculate the thermal conductivity of the mixture in a control volume
double Species_Solver::JANAF_ThermalCond(double T, int i, int j, int k){
int i;
double ThermalCond;
double lambda = 0.0;

    for (i = 0; i < N_Species; i++){
        ThermalCond = exp(Species[SP].lambda_coeff[0] + Species[SP].lambda_coeff[1] * log(T) + Species[SP].lambda_coeff[2] * pow(log(T), 2) + Species[SP].lambda_coeff[3] * pow(log(T), 3));
        lambda += Species[SP].Y_Pres[LM(i,j,k,0)] * ThermalCond;
    }

    return lambda;

}