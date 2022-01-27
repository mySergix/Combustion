//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MATRIX INDEX FUNCTIONS                                //
//------------------------------------------------------------------------------------------------//

// Local Index of Collocated and Staggered Matrix

// Local Index Collocated mesh
#define LP(i,j,k,dim) (NY + 2*HP)*(NZ + 2*HP)*((i) - Ix[Rango] + HP) + ((j) + HP)*(NZ + 2*HP) + ((k) + HP) + ((Fx[Rango] - Ix[Rango] + 2*HP)*(NY + 2*HP)*(NZ + 2*HP) * dim)

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
#define BOTTOM(i,j,k) NZ*((i) - Ix[Rango] + 1) + (k)

// Local Index Top Side
#define TOP(i,j,k) NZ*((i) - Ix[Rango] + 1) + (k)

// Local Index Here Side
#define HERE(i,j,k) NY*((i) - Ix[Rango] + 1) + (j)

// Local Index There Side
#define THERE(i,j,k) NY*((i) - Ix[Rango] + 1) + (j) 


// Global Index of Collocated and Staggered Matrix

// Global Index Collocated mesh
#define GP(i,j,k,Dim) (NY + 2*HP)*(NZ + 2*HP)*((i) + HP)  + ((j) + HP)*(NZ + 2*HP) + ((k) + HP) + (NY + 2*HP)*(NZ + 2*HP)*(NX + 2*HP)*(Dim)

// Global Index Staggered U mesh
#define GU(i,j,k,Dim) (NY + 2*Halo)*(NZ + 2*Halo)*((i) + Halo)  + ((j) + Halo)*(NZ + 2*Halo) + ((k) + Halo)

// Global Index Staggered V mesh
#define GV(i,j,k,Dim) (NY + 2*Halo + 1)*(NZ + 2*Halo)*((i) + Halo)  + ((j) + Halo)*(NZ + 2*Halo) + ((k) + Halo)

// Global Index Staggered W mesh
#define GW(i,j,k,Dim) (NY + 2*Halo)*(NZ + 2*Halo + 1)*((i) + Halo)  + ((j) + Halo)*(NZ + 2*Halo + 1) + ((k) + Halo)
