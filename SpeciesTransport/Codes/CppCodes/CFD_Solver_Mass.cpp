//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR N-S MASS EQUATION CALCULATIONS                        //
//------------------------------------------------------------------------------------------------//

// Function to calculate the step contribution to the density
void CFD_Solver::Get_StepContribution_Density(Property &PropertyName){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                PropertyName.ContributionPres[LM(i,j,k,0)] = - PropertyName.Convective[LM(i,j,k,0)];
            }
        }
    }

}