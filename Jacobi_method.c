#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int diag_domin(double** A, int n)
{
    int i, j;
    double sum;

    for (i = 0; i < n; i++)
    {
        sum = 0.0;
        for (j = 0; j < n; j++)
        {
            if (j != i)
            {
                sum += fabs(A[i][j]);
            }
        }
        if (fabs(A[i][i]) <= sum)
        {
            return 0; // la matrice n'est pas diagonalement dominante
        }
    }
    return 1; // la matrice est diagonalement dominante
}

void jacobi(double** A, double* b, double* x, int n)
{
    int i, j, k, aide;
    double* x_new;
    double diff, sum;

    x_new = (double*)malloc(n * sizeof(double));


    for (k = 1; k <= 1000; k++)
    {
        aide = 1; // indicateur de convergence

        // calculer les nouvelles valeurs de x
        for (i = 0; i < n; i++)
        {
            sum = 0.0;
            for (j = 0; j < n; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        // v�rifier si la solution a converg�
        for (i = 0; i < n; i++)
        {
            diff = fabs(x_new[i] - x[i]);
            if (diff > EPSILON)
            {
                aide = 0;
                break;
            }
        }

        // remplacer x par les nouvelles valeurs
        for (i = 0; i < n; i++)
        {
            printf("sol = %f\n", x_new[i]); 
            x[i] = x_new[i];
        }

        // sortir de la boucle
        if (aide)
        {
            printf("Methode de Jacobi converge en %d iterations.\n", k);
            break;
        }
    }

    free(x_new);
}


int main()
{

    double** A;
    double* b;
    double* x;
    int i, j,n;

    printf("Entrez la taille du systeme : ");
    scanf("%d", &n);

    // allocation de m�moire pour la matrice A, le vecteur b et le vecteur x
    A = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++)
    {
        A[i] = (double*)malloc(n * sizeof(double));
    }
    b = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));

// remplissage
    printf("\n * Remplissage de la matrice A:\n\n");

    for (int i=0;i<n;i++)
    {
        printf(" | ");
        for (int j=0;j<n;j++)
         {
                printf("A[%d][%d] = ",i+1,j+1);
                scanf("%f",&A[i][j]);
        }
    }
     printf("\n * remplissage du vecteur B:\n\n");

    for (int i=0;i<n;i++)
    {
        printf(" | ");
        printf("B[%d] = ",i+1);
        scanf("%f",&b[i]);
    }



    // initialisation de la solution initiale x

    for (i = 0; i < n; i++)
    {
        x[i] = 0.0;
    }

    // v�rification de la diagonale dominante de la matrice A
    if (!diag_domin(A, n))
    {
        printf("La matrice A n'est pas diagonalement dominante.\nOn ne peut pas appliquer la methode de Jacobi\n Essayer une autre methode\n");
        return 1;
    }

    // r�solution du syst�me avec la m�thode de Jacobi
    jacobi(A, b, x, n);

    // affichage de la solution
    printf("La solution est : \n");
    for (i = 0; i < n; i++)
    {
        printf("x[%d] = %f\n", i, x[i]);
    }
    return 0;
}
