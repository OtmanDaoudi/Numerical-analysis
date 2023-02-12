#include<stdio.h>
#include<stdlib.h>

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

void eliminateEntry(float** A, float** B, unsigned dim, unsigned line, unsigned column)
{
    //TODO: function to select a pivot [non-zero]
    float pivot = A[column][column];
    float multiplier = A[line][column]/pivot; 
    for(int col = 0; col < dim; col++) A[line][col] = A[line][col] - multiplier * A[column][col]; //alter A
    B[line][0] = B[line][0] - multiplier * B[column][0]; //alter second memeber   
}

float* gaussianElimination(float** A, float** B, unsigned dim)
{
    //elimination phase
    for(int column = 0; column < dim-1 ; column++)
    {
        for(int line = column+1; line < dim; line++)
        {
            if( *(*(A + line) + column) != 0) eliminateEntry(A, B, dim, line, column); 
        }
    }

    //display current res before final mod
    showMatrix(A, dim, dim); 
    showMatrix(B, dim, 1); 

    //backsolving 
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

int main(void)
{
    unsigned dim; 
    printf("dim(A) = "); 
    scanf("%u", &dim); 

    float* A[dim]; 
    fillMatrix(A, dim, dim); 
    showMatrix(A, dim, dim); 

    float* B[dim];
    fillMatrix(B, dim, 1); 
    showMatrix(B, dim, 1); 

    float* res = gaussianElimination(A, B, dim);
    // showMatrix(&res, 1, dim);  
    printf("x = %f | y = %f | z = %f\n", res[0], res[1], res[2]); 
    return 0; 
}