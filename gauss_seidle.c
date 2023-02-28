#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 100
#define EPSILON 10^(-6)

void gauss_seidel(double A[N][N], double b[N], double x[N], int n)
{
    int i, j, k;
    double norme, somme;
    double y[N];
    for (k = 0; k < 1000; k++)   // 1000 nombre maximal d'it�rations
    {
        norme = 0.0;


        for (i = 0; i < n; i++)
        {
            somme = b[i];
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    somme -= A[i][j] * x[j];
                }
            }
            x[i] = somme / A[i][i];
            norme += fabs(x[i] - y[i]);
        }


        if (norme < EPSILON)
        {
            printf("Gauss-Seidel a converg� en %d it�rations\n", k + 1);
            return;
        }
    }

    printf("Gauss-Seidel n'a pas converg� en %d it�rations\n", k);
}

int main()
{
    double A[N][N], b[N], x[N];
    int n, i, j;

    printf("Entrez la taille de la matrice : ");
    scanf("%d", &n);

    printf("Entrez les coefficients de la matrice A :\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("A[%d][%d]= ",i+1,j+1);
            scanf("%lf", &        A[i][j]);
        }
    }


    printf("Entrez les coefficients du vecteur b :\n");
    for (i = 0; i < n; i++)
    {
        scanf("%lf", &b[i]);
    }




//initier x � 0

    for (i = 0; i < n; i++)
    {
        x[i] = 0.0;
    }

// r�soudre le syst�me avec la m�thode de Gauss-Seidel
    gauss_seidel(A, b, x, n);

//resultat

    printf("Solution :\n");
    for (i = 0; i < n; i++)
    {
        printf("%f\n", x[i]);
    }

    return 0;
}
