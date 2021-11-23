//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR POST PROCESSING CLASS                                 //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcess.h"

using namespace std;

#define LM(i,j,k,dim) ((NY + 2*Halo) * (NZ + 2*Halo)) * ((i) - Ix[Rango] + Halo) + ((NZ + 2*Halo) * ((j) + Halo)) + ((k) + Halo) + ((Fx[Rango] - Ix[Rango] + 2*Halo) * (NY + 2*Halo) * (NZ + 2*Halo)) * (dim)
#define GM(i,j,k,dim) (NY + 2*Halo) * (NZ + 2*Halo) * ((i) + Halo) + (NZ + 2*Halo) * ((j) + Halo) + ((k) + Halo) + (NY + 2*Halo) * (NZ + 2*Halo) * (NX + 2*Halo) * (dim)

PostProcess::PostProcess(Memory M1, ReadData R1, Parallel P1){
	
    // Data necessary

        // Meshing data
        NX = R1.ProblemNumericalData[0];
        NY = R1.ProblemNumericalData[1];
        NZ = R1.ProblemNumericalData[2];

        // Parallel computing data
        Rango = P1.Rango;
        Procesos = P1.Procesos;

        Halo = 2;

        Ix = M1.AllocateInt(Procesos, 1, 1, 1);
        Fx = M1.AllocateInt(Procesos, 1, 1, 1);

        for (int i = 0; i < Procesos; i++){
            Ix[i] = P1.Ix[i];
            Fx[i] = P1.Fx[i];
        }

}

//Pasar los resultados de un Scalar a un VTK en 3D
void PostProcess::GlobalEscalarVTK(Mesher MESH, string Carpeta, string Variable, string NombreFile, double *ScalarMatrix, int HALO){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<"../ParaviewResults/"<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<Variable<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<(NX + 2*HALO)<<"   "<<(NY + 2*HALO)<<"   "<<(NZ + 2*HALO)<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<(NX + 2*HALO) * (NY + 2*HALO) * (NZ + 2*HALO)<<"   "<<"double"<<endl;
	
	for(k = - HALO; k < NZ + HALO; k++){
		for(j = - HALO; j < NY + HALO; j++){
			for(i = - HALO; i < NX + HALO; i++){
				file<<MESH.Global_Node_Mesh[GM(i,j,k,0)]<<"   "<<MESH.Global_Node_Mesh[GM(i,j,k,1)]<<"   "<<MESH.Global_Node_Mesh[GM(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
	file<<"POINT_DATA"<<"   "<<(NX + 2*HALO) * (NY + 2*HALO) * (NZ + 2*HALO)<<endl;
    file<<"SCALARS "<<Variable<<" double"<<endl;
    file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
    file<<endl;
	for(k = - HALO; k < NZ + HALO; k++){
		for(j = - HALO; j < NY + HALO; j++){
			for(i = - HALO; i < NX + HALO; i++){
				file<<ScalarMatrix[GM(i,j,k,0)]<<" ";
			}
		}
	}

    file.close();

}

//Pasar los resultados de un Scalar a un VTK en 3D
void PostProcess::GlobalVectorialVTK(Mesher MESH, string Carpeta, string Variable, string NombreFile, Global &GlobalMatrix, int HALO){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<"../ParaviewResults/"<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<Variable<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<(NX + 2*HALO)<<"   "<<(NY + 2*HALO)<<"   "<<(NZ + 2*HALO)<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<(NX + 2*HALO) * (NY + 2*HALO) * (NZ + 2*HALO)<<"   "<<"double"<<endl;
	
	for(k = - HALO; k < NZ + HALO; k++){
		for(i = - HALO; i < NX + HALO; i++){
			for(j = - HALO; j < NY + HALO; j++){
				file<<MESH.Node_Mesh[GM(i,j,k,0)]<<"   "<<MESH.Node_Mesh[GM(i,j,k,1)]<<"   "<<MESH.Node_Mesh[GM(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
    file<<"POINT_DATA"<<"   "<<(NX + 2*HALO) * (NY + 2*HALO) * (NZ + 2*HALO)<<endl;
    file<<"VECTORS "<<Variable<<" double"<<endl;
    file<<endl;

	for(k = - HALO; k < NZ + HALO; k++){
		for(i = - HALO; i < NX + HALO; i++){
			for(j = - HALO; j < NY + HALO; j++){	
				file<<GlobalMatrix.U[GM(i,j,k,0)]<<" "<<GlobalMatrix.V[GM(i,j,k,0)]<<" "<<GlobalMatrix.W[GM(i,j,k,0)]<<endl;
			}
		}
	}

    file.close();

}