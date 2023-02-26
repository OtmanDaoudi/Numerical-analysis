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

//cholseky + crammer

// int main(){
//     /// Test de Crammer
//     float d;
//     float **M,**L,*Y,*X,B[n]={0.5,0.25,1},A[n][n]={ {1,0,-0.25,-0.25},
//                         {0,1,-0.25,-0.25},
//                         {-0.25,-0.25,1,-0.25}}; ///{1,0,-0.25},{0,1,-0.25},{-0.25,-0.25,1}
//     M = (float**)calloc(n,sizeof(float*));
//     for(int i=0;i<n;i++)
//      {
//             M[i] = (float *)calloc(n,sizeof(float));
//      }
//     copyMat(A,M,n); //A=M
//     afficher_matStat(A,n);
//     //d=determinant(A[m],n);
//     //afficher_matStat(A,n); ///appel erron� de la fct
//     //printf("det(A)=%f",d);
//     //d=determinantV2(A,n); ///M triangulaire par precess de gauss
//     //afficher_matStat(A,n);
//     //afficher_mat(M,n);
//     //printf("det(A)=%f",d);
//     //afficher_mat(M,n);
//     //echangeCol(M,B,1,n);
//     printf("\naffichage du second membre :\n");
//     afficher_vect(B,n);
//     printf("\nsolution par methode de Crammer:\n");
//     afficher_unk(Crammer(A,B,n),n);
//     free(M);

//     ///-------------------------------------------------------------------------------
//     /// Test de Cholesky


//     /// Cholesky Manuelle:
//     /* L = cholesky(A,n);
//     afficher_matStat(A,n);
//     afficher_mat(L,n); ///L Triangulaire inferieur
//     afficher_mat(transposeDyn(L,n)); ///Transpose(L) triangulaire sup
//     Y = reso_sysTriDyn(transposeDyn(L,n),B,n);
//     afficher_unk(Y,n);
//     X = reso_sysTriDynInv(L,Y,n);
//     afficher_unk(X,n); */
//     ///Cholesky Auto:
//     printf("\nsolution par methode de cholesky:\n");
//     afficher_unk(solChol(A,B,n),n);
//     return 0;
//     return 0;
// }
