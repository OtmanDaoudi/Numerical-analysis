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

void eliminateEntry(float** A, unsigned dim, unsigned line, unsigned column, float** L)
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
    L[line][column] = multiplier; 
    // B[line][0] = B[line][0] - multiplier * B[column][0]; //alter second memeber   
}


float** LU(float** A,unsigned dim) //return L
{
    //initialise L
    float** L = (float**)malloc(sizeof(float*) * dim); 
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
            if( *(*(A + line) + column) != 0) eliminateEntry(A, dim, line, column, L); 
        }
    }
    //A = U !!
    return L;    
}

//LU decomposition + solution
int main(void)
{
    unsigned dim; 
    printf("dim(A) = "); 
    scanf("%u", &dim); 

    float* A[dim]; 
    fillMatrix(A, dim, dim); 
    printf("A : \n"); 
    showMatrix(A, dim, dim); 

    float** L = LU(A, dim); 

    printf("U : \n"); 
    showMatrix(A, dim, dim);

    printf("L : \n"); 
    showMatrix(L, dim, dim); 


    return 0; 
}