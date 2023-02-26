#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define n 3

//gaussian elimination
void    fillMatrix(float** matrix, unsigned lines, unsigned columns); 
void    showMatrix(float** matrix, unsigned lines, unsigned columns); 
float   eliminateEntry(float** A, unsigned dim, unsigned line, unsigned column); 
float*  gaussianElimination(float** A, float** B, unsigned dim); 
float*  backSolve(float** A, float** B, unsigned dim); 
float*  backSolveReversed(float** A, float** B, unsigned dim); //for Ls 

//LU decomposition
float* LU(float** A,float** B, unsigned dim, float** L); 

//Crammer + Cholesky
float   determinant(float A[][n],float **mat_res, int dim); 
float   determinantV2(float A[][n], int dim); 
void    copyMat(float A[n][n],float **mat_cpy,int dim); 
void    echangeCol(float **A,float *B,int col,int dim); 
float*  Crammer(float **A,float *B,int dim); 
float** cholesky(float A[n][n],int dim);
float** transposeDyn(float **L,int dim);
float*  reso_sysTriDyn(float **L,float B[n],int dim);
float*  reso_sysTriDynInv(float **L,float *B,int dim);
float*  solChol(float A[][n],float B[],int dim);
void    afficher_vect(float *B,int dim);
void    afficher_unk(float *X,int dim);
void    afficher_mat(float **X,int dim);
void    afficher_matStat(float X[n][n],int dim);