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


// //LU main function
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

//     float* L[dim];
//     float* solution = LU(A, B, dim, L); 

//     printf("U : \n"); 
//     showMatrix(A, dim, dim);

//     printf("L : \n"); 
//     showMatrix(L, dim, dim); 

//     printf("Solution : \n");
//     showMatrix(&solution, 1, dim);

//     return 0; 
// }

//cholseky + crammer

// int main(){
//     /// Test de Crammer
//     float **M,B[n]={0.5,0.25,1},A[n][n]={ {1,0,-0.25},
//                         {0,1,-0.25},
//                         {-0.25,-0.25,1}}; ///{1,0,-0.25},{0,1,-0.25},{-0.25,-0.25,1}
//     for(int i=0;i<n;i++)
//      {
//             M[i] = (float *)calloc(n,sizeof(float));
//      }
//      //A=M
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
//     afficher_unk(Crammer((float**)A,B,n),n);
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


int main(void)
{
    unsigned option; 
    do
    {
        printf("========= Menu ==========\n");
        printf("Methodes directes: \n"); 
        printf("\t1-Method LU\n\t2-Method de gauss\n\t3-Methode de crammer\n\t4-Methode LLt\n\n");  
        printf("Methodes iteratives: \n"); 
        printf("\t5-Methode de jacobi\n\t6-Methode de Gausse-Sciedle\n");  
        printf("7-Quitter"); 
        printf("=========================\n");
        printf("Option : "); scanf("%u", &option); 

        unsigned dim; 
        printf("dim(A) = "); scanf("%u", &dim); 
        
        float *A[dim], *B[dim], AAux[dim][dim], Baux[dim];
        // M = (float**)calloc(dim,sizeof(float*));

        printf("Matrice A : \n"); 
        fillMatrix(A, dim, dim); 
        printf("\nMatrice B : \n"); 
        fillMatrix(B, dim, 1);
    
        float* L[dim], *solution;

        switch(option)
        {
            case 3: //Methode de crammer
                for(int i=0; i<dim; i++) Baux[i] = B[i][0]; 
                for(int l=0; l<dim; l++) for(int c=0; c<dim; c++) AAux[l][c] = A[l][c]; 
                printf("\nSolution : \n"); 
                afficher_unk(Crammer(AAux,Baux,dim),dim);  
                break;
            case 4: //Methode LLt
                for(int i=0; i<dim; i++) Baux[i] = B[i][0]; 
                for(int l=0; l<dim; l++) for(int c=0; c<dim; c++) AAux[l][c] = A[l][c]; 
                printf("\nSolution : \n"); 
                afficher_unk(solChol(AAux,Baux,dim),dim);
                break; 
            case 1: //Methode LU
                solution = LU(A, B, dim, L);
                printf("\nU : \n"); 
                showMatrix(A, dim, dim);

                printf("L : \n"); 
                showMatrix(L, dim, dim); 

                printf("Solution : \n");
                showMatrix(&solution, 1, dim);
                break; 
            case 2: //Methode de gausse
                solution = gaussianElimination(A, B, dim);
                printf("\nMatrice A apres les operations : \n"); 
                showMatrix(A, dim, dim); 

                printf("Matrice B apres les operations : \n"); 
                showMatrix(B, dim, 1); 

                printf("\nSolution : \n"); 
                showMatrix(&solution, 1, dim);  
                break; 
            case 5: //methode jacobi
                for(int i=0; i<dim; i++) Baux[i] = B[i][0]; 
                solution = (float*)malloc(sizeof(float) * dim);
                printf("Solution intitiale : \n"); 
                for(int i=0; i<dim; i++)
                {
                    printf("Mat[%d] = ", i); 
                    scanf("%f", solution+i);
                }
                jacobi(A, Baux, solution, dim);
                printf("\nSolution : \n"); 
                showMatrix(&solution, 1, dim);
                break; 
            case 6: //gauss siedle
                for(int i=0; i<dim; i++) Baux[i] = B[i][0]; 
                solution = (float*)malloc(sizeof(float) * dim);
                printf("Solution intitiale : \n"); 
                for(int i=0; i<dim; i++)
                {
                    printf("Mat[%d] = ", i); 
                    scanf("%f", solution+i);
                }
                gauss_seidel(A, Baux, solution, dim); 
                printf("\nSolution : \n"); 
                showMatrix(&solution, 1, dim);
                break; 
        }
    } while (option != 7);
}