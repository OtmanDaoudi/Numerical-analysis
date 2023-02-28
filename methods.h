#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define EPSILON 1e-6

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
float   determinant(float** A,float **mat_res, int dim); 
float   determinantV2(float** A, int dim); 
void    copyMat(float** A,float **mat_cpy,int dim); 
void    echangeCol(float **A,float *B,int col,int dim); 
float*  Crammer(float **A,float *B,int dim); 
float** cholesky(float** A,int dim);
float** transposeDyn(float **L,int dim);
float*  reso_sysTriDyn(float **L,float* B,int dim);
float*  reso_sysTriDynInv(float **L,float *B,int dim);
float*  solChol(float** A,float B[],int dim);
void    afficher_vect(float *B,int dim);
void    afficher_unk(float *X,int dim);
void    afficher_mat(float **X,int dim);
void    afficher_matStat(float** X,int dim);

//jacobi + gauss sciedle
int diag_domin(float** A, int n); 
void jacobi(float** A, float* b, float* x, int n); 
void gauss_seidel(double** A, double* b, double* x, int n); 
