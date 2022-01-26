//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MATRIX INDEX FUNCTIONS                                //
//------------------------------------------------------------------------------------------------//

// Local Index of Collocated and Staggered Matrix

// Local Index Collocated mesh
#define LP(i,j,k,dim) (NY + 2*HP)*(NZ + 2*HP)*((i) - Ix[Rango] + HP) + (j + HP)*(NZ + 2*HP) + ((k) + HP) + ((Fx[Rango] - Ix[Rango] + 2*HP)*(NY + 2*HP)*(NZ + 2*HP) * dim)

// Local Index Staggered U mesh
#define LU(i,j,k,dim) (NY + 2*Halo)*(NZ + 2*Halo)*((i) - Ix[Rango] + Halo) + (NZ + 2*Halo)*((j) + Halo) + ((k) + Halo) + ((Fx[Rango] - Ix[Rango] + 2*HP + 1)*(NY + 2*HP)*(NZ + 2*HP) * dim)

// Local Index Staggered V mesh
#define LV(i,j,k,dim) (NY + 2*Halo + 1)*(NZ + 2*Halo)*((i) - Ix[Rango] + Halo) + (NZ + 2*Halo)*((j) + Halo) + ((k) + Halo) +  ((Fx[Rango] - Ix[Rango] + 2*HP)*(NY + 2*HP + 1)*(NZ + 2*HP) * dim)

// Local Index Staggered W mesh
#define LW(i,j,k,dim) (NY + 2*Halo)*(NZ + 2*Halo + 1)*((i) - Ix[Rango] + Halo) + (NZ + 2*Halo + 1)*((j) + Halo) + ((k) + Halo) + ((Fx[Rango] - Ix[Rango] + 2*HP)*(NY + 2*HP)*(NZ + 2*HP + 1) * dim)   

// Local Index Poisson Coefficients
#define LA(i,j,k,dim) (NY * NZ)*((i) - Ix[Rango]) + (NZ)*(j) + (k)


// Boundary Conditions Index

// Local Index Left Side
#define LEFT(i,j,k) (NZ*(j)) + (k)

// Local Index Right Side
#define RIGHT(i,j,k) (NZ*(j)) + (k)

// Local Index Bottom Side
#define BOTTOM(i,j,k) NZ*((i) - Ix[Rango]) + (k)

// Local Index Top Side
#define TOP(i,j,k) NZ*((i) - Ix[Rango]) + (k)

// Local Index Here Side
#define HERE(i,j,k) NY*((i) - Ix[Rango]) + (j)

// Local Index There Side
#define THERE(i,j,k) NY*((i) - Ix[Rango]) + (j) 


// -----------------------------------------------------------------------------------------------------------------------------------------------------------
//Global Index P Mesh
#define GP(i,j,k,Dim) (NY + 2*HaloPressure)*(NZ + 2*HaloPressure)*((i) + HaloPressure) + ((j) + HaloPressure) + ((k) + HaloPressure)*(NY + 2*HaloPressure) + (NX + 2*HaloPressure)*(NY + 2*HaloPressure)*(NZ + 2*HaloPressure)*(Dim)

//Global Index U Mesh
#define GU(i,j,k,Dim) (NY + 2*HaloU)*(NZ + 2*HaloU)*((i) + HaloU) + ((j) + HaloU) + ((k) + HaloU)*(NY + 2*HaloU) + (NX + 1 + 2*HaloU)*(NY + 2*HaloU)*(NZ + 2*HaloU)*(Dim)

//Global Index V Mesh
#define GV(i,j,k,Dim) (NY + 1 + 2*HaloV)*(NZ + 2*HaloV)*((i) + HaloV) + ((j) + HaloV) + ((k) + HaloV)*(NY + 1 + 2*HaloV) + (NX + 2*HaloV)*(NY + 1 + 2*HaloV)*(NZ + 2*HaloV)*(Dim)

//Global Index W Mesh
#define GW(i,j,k,Dim) (NY + 2*HaloW)*(NZ + 1 + 2*HaloW)*((i) + HaloW) + ((j) + HaloW) + ((k) + HaloW)*(NY + 2*HaloW) + (NX + 2*HaloW)*(NY + 2*HaloW)*(NZ + 1 + 2*HaloW)*(Dim)
