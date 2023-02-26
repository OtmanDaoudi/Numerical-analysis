#include "methods.h"

// //gaussian elimination main function

// int main(void)
// {
//     unsigned dim; 
//     printf("dim(A) = "); 
//     scanf("%u", &dim); 

//     float* A[dim]; 
//     printf("A : \n"); 
//     fillMatrix(A, dim, dim); 
//     showMatrix(A, dim, dim); 

//     float* B[dim];
//     printf("B : \n"); 
//     fillMatrix(B, dim, 1); 
//     showMatrix(B, dim, 1); 

//     float* res = gaussianElimination(A, B, dim);
//     showMatrix(&res, 1, dim);  
//     // printf("x = %f | y = %f | z = %f\n", res[0], res[1], res[2]); 
//     return 0; 
// }


//LU main function
int main(void)
{
    unsigned dim; 
    printf("dim(A) = "); 
    scanf("%u", &dim); 

    float* A[dim];
    printf("A : \n");  
    fillMatrix(A, dim, dim); 
    showMatrix(A, dim, dim); 

    float* B[dim];
    printf("B : \n");  
    fillMatrix(B, dim, 1); 
    showMatrix(B, dim, 1);

    float* L[dim];
    float* solution = LU(A, B, dim, L); 

    printf("U : \n"); 
    showMatrix(A, dim, dim);

    printf("L : \n"); 
    showMatrix(L, dim, dim); 

    printf("Solution : \n");
    showMatrix(&solution, 1, dim);

    return 0; 
}