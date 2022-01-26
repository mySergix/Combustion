//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MESHER MEMORY FUNCTIONS                               //
//------------------------------------------------------------------------------------------------//

// Memory allocation for each matrix
void Mesher::Allocate_MesherMemory(Memory M1){

    // Code
	// 0 -> Coordenada X
	// 1 -> Coordenada Y
	// 2 -> Coordenada Z

	// Meshes of the simulation
	MP = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 3); // Collocated mesh
	MU = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo + 1, NY + 2*Halo, NZ + 2*Halo, 3); // Staggered U mesh
	MV = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo + 1, NZ + 2*Halo, 3); // Staggered V mesh
	MW = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, NZ + 2*Halo + 1, 3); // Staggered W mesh

    // Meshes nodal distances 
	DeltasMP = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 3); // Collocated mesh distances
	DeltasMU = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo + 1, NY + 2*Halo, NZ + 2*Halo, 3); // Staggered U mesh distances
	DeltasMV = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo + 1, NZ + 2*Halo, 3); // Staggered V mesh distances
	DeltasMW = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, NZ + 2*Halo + 1, 3); // Staggered W mesh distances
	
    // Meshes surfaces of the CV
	SupMP = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 3); // Collocated mesh surfaces
	SupMU = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo + 1, NY + 2*Halo, NZ + 2*Halo, 3); // Staggered U mesh surfaces
	SupMV = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo + 1, NZ + 2*Halo, 3); // Staggered V mesh surfaces
	SupMW = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, NZ + 2*Halo + 1, 3); // Staggered W mesh surfaces

    // Volumes of the CV of the meshes
	VolMP = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 3); // Collocated mesh CV Volumes
	VolMU = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo + 1, NY + 2*Halo, NZ + 2*Halo, 3); // Staggered U mesh CV Volumes
	VolMV = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo + 1, NZ + 2*Halo, 3); // Staggered V mesh CV Volumes
	VolMW = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*Halo, NY + 2*Halo, NZ + 2*Halo + 1, 3); // Staggered W mesh CV Volumes

}

// Function to delete all the memory from the mesher class
void Mesher::Delete_MesherMemory(){

	// Meshes of the simulation
	delete [] MP; // Collocated mesh
	delete [] MU; // Staggered U mesh
	delete [] MV; // Staggered V mesh
	delete [] MW; // Staggered W mesh

    // Meshes nodal distances 
	delete [] DeltasMP; // Collocated mesh distances
	delete [] DeltasMU; // Staggered U mesh distances
	delete [] DeltasMV; // Staggered V mesh distances
	delete [] DeltasMW; // Staggered W mesh distances
	
    // Meshes surfaces of the CV
	delete [] SupMP; // Collocated mesh surfaces
	delete [] SupMU; // Staggered U mesh surfaces
	delete [] SupMV; // Staggered V mesh surfaces
	delete [] SupMW; // Staggered W mesh surfaces

    // Volumes of the CV of the meshes
	delete [] VolMP; // Collocated mesh CV Volumes
	delete [] VolMU; // Staggered U mesh CV Volumes
	delete [] VolMV; // Staggered V mesh CV Volumes
	delete [] VolMW; // Staggered W mesh CV Volumes

	delete [] Ix;
	delete [] Fx;
	
}