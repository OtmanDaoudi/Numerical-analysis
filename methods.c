#include "methods.h"

void fillMatrix(float** matrix, unsigned lines, unsigned columns)
{
    if(lines == 0 || columns == 0) return; 
    for(int line=0; line<lines; line++)
    {
        *(matrix + line) = (float*)malloc(sizeof(float) * columns); 
        for(int column=0; column<columns; column++)
        {
            printf("Mat[%d][%d] = ", line, column); 
            scanf("%f", *(matrix + line) + column); 
        }
    }
}

void showMatrix(float** matrix, unsigned lines, unsigned columns)
{
    if(matrix == NULL || *matrix == NULL || lines == 0 || columns == 0) return; 
    printf("\n"); 
    for(int line=0; line<lines; line++)
    {
        for(int column=0; column<columns; column++) printf("%f\t", *(*(matrix + line) + column)); 
        printf("\n"); 
    }    
    printf("\n"); 
}

float eliminateEntry(float** A, unsigned dim, unsigned line, unsigned column)
{
    //TODO: function to select a pivot [non-zero]
    float pivot = A[column][column];
    if(pivot == 0) 
    {
        printf("pivot is 0\n"); 
        exit(EXIT_FAILURE); 
    }
    float multiplier = A[line][column]/pivot; 
    for(int col = 0; col < dim; col++) A[line][col] = A[line][col] - multiplier * A[column][col]; //alter A
    return multiplier; 
}

float* backSolve(float** A, float** B, unsigned dim)
{
    float* resVector = (float*)malloc(sizeof(float) * dim); 
    for(int line = dim - 1;  line >= 0 ; line--)
    {
        resVector[line] = B[line][0]; 
        for(int column = line+1; column <= dim - 1 ; column++)
        {
            resVector[line] += - A[line][column] * resVector[column];
        }
        resVector[line] /= A[line][line];
    }
    return resVector; 
}

float* backSolveReversed(float** A, float** B, unsigned dim)
{
    float* resVector = (float*)malloc(sizeof(float) * dim); 
    for(int line = 0;  line < dim ; line++)
    {
        resVector[line] = B[line][0]; 
        for(int column = line - 1; column >= 0 ; column--)
        {
            resVector[line] += - A[line][column] * resVector[column];
        }
        resVector[line] /= A[line][line];
    }
    return resVector;
}

float* gaussianElimination(float** A, float** B, unsigned dim)
{
    //elimination phase
    for(int column = 0; column < dim-1 ; column++)
    {
        for(int line = column+1; line < dim; line++)
        {
            if( *(*(A + line) + column) != 0) 
            {
                B[line][0] = B[line][0] - eliminateEntry(A, dim, line, column) * B[column][0]; //alter second memeber   
            }
        }
    }

    //TODO: remove showMats
    //display current res before final mod
    showMatrix(A, dim, dim); 
    showMatrix(B, dim, 1); 

    //backsolving 
    return backSolve(A, B, dim);    
}

float* LU(float** A,float** B, unsigned dim, float** L) //return L
{
    //initialise L
    // float** L = (float**)malloc(sizeof(float*) * dim); 
    for(int line = 0; line < dim ; line++)
    {
        L[line] = (float*)malloc(sizeof(float) * dim); 
        for(int column = 0; column < dim ; column++)
        {
            if(line == column) L[line][column] = 1; //diagonal 
            else if(line < column) L[line][column] = 0; //L has 0s above the diagonal
        }
    }
    //elimination phase
    for(int column = 0; column < dim-1 ; column++)
    {
        for(int line = column+1; line < dim; line++)
        {
            if( *(*(A + line) + column) != 0) L[line][column] = eliminateEntry(A, dim, line, column); 
        }
    }
    //A = U !!
    //backsolving phase
    //Ax=B ==> LUx = B
    /*
        backsolving on two stges:
            Lc=B
            Ux=c
    */

    float* c = backSolveReversed(L, B, dim); 
    float* tempC[dim]; 
    for(int line=0; line<dim; line++)
    {
        tempC[line] = (float*)malloc(sizeof(float)); 
        tempC[line][0] = c[line]; 
    }
    
    return backSolve(A, tempC, dim); 
    return NULL; 
}
