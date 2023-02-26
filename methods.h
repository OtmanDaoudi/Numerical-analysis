#include<stdio.h>
#include<stdlib.h>

//gaussian elimination
void fillMatrix(float** matrix, unsigned lines, unsigned columns); 
void showMatrix(float** matrix, unsigned lines, unsigned columns); 
float eliminateEntry(float** A, unsigned dim, unsigned line, unsigned column); 
float* gaussianElimination(float** A, float** B, unsigned dim); 
float* backSolve(float** A, float** B, unsigned dim); 
float* backSolveReversed(float** A, float** B, unsigned dim); //for Ls 

//LU decomposition
float* LU(float** A,float** B, unsigned dim, float** L); 
