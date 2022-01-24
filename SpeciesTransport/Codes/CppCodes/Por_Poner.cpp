 {

        

        

        // Velocities
        //Get_Pressure();
        //Get_PressureGradient(MESH);

        // Velocity U
        Get_ConvectiveTerm(MESH, U);
        Get_StepContribution_Velocity(U, Pressure.Gradient_X);

        // Velocity V
        Get_ConvectiveTerm(MESH, V);
        Get_StepContribution_Velocity(V, Pressure.Gradient_Y);

        // Velocity W
        Get_ConvectiveTerm(MESH, W);
        Get_StepContribution_Velocity(W, Pressure.Gradient_Z);

        Get_TemporalIntegration(U);
        Get_TemporalIntegration(V);
        Get_TemporalIntegration(W);

        if (Step % 100 == 0){
            Get_ConvergenceCriteria();
            P1.SendMatrixToZero(Density.Pres, GlobalMatrix.Density);
            
            P1.SendMatrixToZero(U.Pres, GlobalMatrix.U);
            P1.SendMatrixToZero(V.Pres, GlobalMatrix.V);
            P1.SendMatrixToZero(W.Pres, GlobalMatrix.W);

            if (Rango == 0){
                cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiff<<endl;
                sprintf(FileName_1, "Density_Field_Step_%d", Step);
                PP1.GlobalEscalarVTK(MESH, "CombustionResults/", "Densidad", FileName_1, GlobalMatrix.Density, 0);

                sprintf(FileName_1, "Velocities_Field_Step_%d", Step);
                PP1.GlobalVectorialVTK(MESH, "CombustionResults/", "Velocities", FileName_1, GlobalMatrix.U, GlobalMatrix.V, GlobalMatrix.W, 0);
            }

            MPI_Barrier(MPI_COMM_WORLD);	
        }
        
        UpdateField(Density);
        UpdateField(U);
        UpdateField(V);
        UpdateField(W);

    }
    