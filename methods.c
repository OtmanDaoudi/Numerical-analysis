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
}


/// partie Crammer
float determinant(float** A,float **mat_res, int dim){ ///triangularise la mat_res et retourne det(A)
    int signe = 1;
    float pivot, coef,det=1;

   for(unsigned int i = 0; i < dim; i++)
   {
       pivot = mat_res[i][i];
       if(pivot == 0)
       {
           int indice = -1;
           for(unsigned int j = i + 1; j < dim; j++)
           {
               if(mat_res[j][i] != 0)
               {
                   indice = j;
                   break;
               }
           }
           if(indice == -1)
           {
               return 0; /// si le pivot est nul et tous en dessous nuls, le d�terminant = 0;
           }
           for(unsigned int k = i; k < dim; k++)
           {
               float temp = mat_res[i][k];
               mat_res[i][k] = mat_res[indice][k];
               mat_res[indice][k] = temp;
           }
           signe *= -1;
           pivot = mat_res[i][i];
       }//printf("p=%f\n",pivot);
       for(unsigned int j = i + 1; j < dim; j++) ///triangularisation de matrice
       {
           coef = mat_res[j][i] / pivot;
           for(unsigned int k = i; k < dim; k++)
           {
               mat_res[j][k] = mat_res[j][k] - coef * mat_res[i][k];
           }
       }
       det *= pivot;
   }

   return signe * det;
}
float determinantV2(float** A, int dim){   /// version sans matrice result
    float **mat_res=(float**)calloc(dim,sizeof(float*));
    for(int i=0;i<dim;i++)
     {
            mat_res[i] = (float *)calloc(dim,sizeof(float));
     }
     copyMat(A,mat_res,dim);
    int signe = 1;
    float pivot, coef,det=1;

   for(unsigned int i = 0; i < dim; i++)
   {
       pivot = mat_res[i][i];
       if(pivot == 0)
       {
           int indice = -1;
           for(unsigned int j = i + 1; j < dim; j++)
           {
               if(mat_res[j][i] != 0)
               {
                   indice = j;
                   break;
               }
           }
           if(indice == -1)
           {
               return 0; /// si le pivot est nul et tous en dessous nuls, le d�terminant = 0;
           }
           for(unsigned int k = i; k < dim; k++)
           {
               float temp = mat_res[i][k];
               mat_res[i][k] = mat_res[indice][k];
               mat_res[indice][k] = temp;
           }
           signe *= -1;
           pivot = mat_res[i][i];
       }//printf("p=%f\n",pivot);
       for(unsigned int j = i + 1; j < dim; j++) ///triangularisation de matrice
       {
           coef = mat_res[j][i] / pivot;
           for(unsigned int k = i; k < dim; k++)
           {
               mat_res[j][k] = mat_res[j][k] - coef * mat_res[i][k];
           }
       }///afficher_matStat(A,n);
       det *= pivot;
   }
    free(mat_res);
   return signe * det;
}

void copyMat(float** A,float **mat_cpy,int dim){ ///n'oublie pas allocation dynam de mat_cpy
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
            *(*(mat_cpy+i)+j)= A[i][j];
        }
    }
}

void echangeCol(float **A,float *B,int col,int dim){
    for (unsigned int i=0;i<dim;i++){
        *(*(A+i)+col) = B[i];
    }
}

float *Crammer(float **A,float *B,int dim){///resolution de AX=B
    float detA=determinantV2(A,dim),dettemp;
    float *X=(float*)calloc(dim,sizeof(float)); ///creation d'un vecteur nul
    float **temp=(float**)calloc(dim,sizeof(float*)); ///matrice dont les colonnes vont chang�es
    for(int i=0;i<dim;i++){
            temp[i] = (float *)calloc(dim,sizeof(float)); ///chaque ligne i contient le nbr "dim" d'elements
     }
     copyMat(A,temp,dim);
     for(unsigned int i=0;i<dim;i++){
         echangeCol(temp,B,i,dim);
         dettemp = determinant(A,temp,dim);
         //afficher_mat(temp,dim);
         //printf("dettemp=%f\tdetA=%f\n",dettemp,detA);
         X[i] = dettemp/detA;
         copyMat(A,temp,dim);
     }
     free(temp);
     return X;
}

/// Partie de cholesky
float ** cholesky(float** A,int dim){///retourne la matrice L triang inf
    float **L;
    float sum=0;
    ///inisialisation de L
    L = (float **)calloc(dim,sizeof(float*));

    for(int i=0;i<dim;i++)
     {
            *(L+i) = (float *)calloc(dim,sizeof(float));
     }

    **L = 1;
    /// first column
    for(unsigned int k=1;k<dim;k++)
     {
            *(*(L+k))= **(A+k) / **L; ///L[0][k]= A[k][0]/L[0][0];
     }
    for(unsigned int i=1;i<dim;i++)
    {

        for(unsigned int k=0;k<=i-1;k++)
            sum+= pow((*(*(L+i)+k)),2);
        ///les elements diagonaux
        *(*(L+i)+i) = sqrt(*(*(A+i)+i)- sum);
        sum = 0;
        /// line i column j
        for(unsigned int j=i+1;j<dim;j++)///  k (1-->i-1) j(i+1-->n)
        {   for(unsigned int k=0;k<=i-1;k++)
                    sum+= *(*(L+j)+k) * (*(*(L+i)+k));

             *(*(L+j)+i)= (*(*(A+j)+i)-sum)/ *(*(L+i)+i) ; ///L[j][i] = (A[j][i]-sum(L[j][k]*L[i][k] ) )/L[i][i]
        }
        sum=0;
    }
    return L;
}

float **transposeDyn(float **L,int dim){
    float **Lt=(float **)calloc(dim,sizeof(float *));
    for(int i=0;i<dim;i++)
     {
            Lt[i] = (float *)calloc(dim,sizeof(float));
     }
     for(unsigned int i=0;i<dim;i++){
        for(unsigned int j=0;j<dim;j++){
            Lt[i][j] = *(*(L+j)+i);
        }
     }
     return Lt;
}

float* reso_sysTriDyn(float **L,float* B,int dim){ ///resoudre syst triang sup par remont�e
    float *X =(float *)calloc(dim,sizeof(float));
    X[dim-1]=B[dim-1]/ *(*(L+dim-1)+dim-1); ///first index = 0 & last n-1 ! A[n-1][n-1]
    for(int i=dim-1;i>=0;i--)
    {
        float sum=0;
        for (int j=dim-1;j>i;j--)
        {
            sum += *(*(L+i)+j)*X[j]; ///A[i][j]
            X[i] = (B[i]-sum)/ *(*(L+i)+i) ;///A[i][i]
        }
    }
    return X;
}

float* reso_sysTriDynInv(float **L,float *B,int dim){///resoudre syst triang inf par descente
    float *X =(float *)calloc(dim,sizeof(float));
    X[0]=B[0]/ **L; ///first index = 0 & last n-1 ! A[n-1][n-1]
    for(int i=1;i<dim;i++)
    {
        float sum=0;
        for (int j=0;j<i;j++)
        {
            sum += *(*(L+i)+j)*X[j]; ///A[i][j]
            X[i] = (B[i]-sum)/ *(*(L+i)+i) ;///A[i][i]
        }
    }
    return X;
}

float *solChol(float** A,float* B,int dim){///retourne uniquement la solution X
        float **L=(float**)calloc(dim,sizeof(float *)); ///se libere apres le retour de la fct
        float **Lt=(float**)calloc(dim,sizeof(float *));///se libere apres le retour de la fct
        float *Y=(float*)calloc(dim,sizeof(float)); ///inconnue temp
        float *X=(float*)calloc(dim,sizeof(float)); ///incon return�e
        for(int i=0;i<dim;i++)
        {
            *(L+i) = (float *)calloc(dim,sizeof(float));
             Lt[i] = (float *)calloc(dim,sizeof(float));
        }
        L = cholesky(A,dim);
        printf("L :\n"); 
        afficher_mat(L, dim); 
        Lt= transposeDyn(L,dim);
        printf("Lt :\n"); 
        afficher_mat(Lt, dim); 
        Y = reso_sysTriDyn(transposeDyn(L,dim),B,dim);
        X = reso_sysTriDynInv(L,Y,dim);
        return X;
}

/// Partie en commun

void afficher_vect(float *B,int dim){ ///afficher vect B
    printf("B=");
    for(int i=0;i<dim;i++)
    {
        if(i==0)
            printf("[%f,",B[i]);
        else if(i == dim-1)
            printf("%f]\n",B[i]);
        else
            printf("%f,",B[i]);
    }
}

void afficher_unk(float *X,int dim){ ///afficher vect X
    printf("X=");
    for(int i=0;i<dim;i++)
    {
        if(i==0)
            printf("[%f,",X[i]);
        else if(i == dim-1)
            printf("%f]\n",X[i]);
        else
            printf("%f,",X[i]);
    }
}

void afficher_mat(float **X,int dim){///affichage de matrice dynamique
    // printf("L=\n");
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
        if(j==0)
            printf("[%f,",*(*(X+i)+j));
        else if(j == dim-1)
            printf("%f]\n",*(*(X+i)+j));
        else
            printf("%f,",*(*(X+i)+j));
        }
    }
}

void afficher_matStat(float** X,int dim){///affichage de matrice statique
    printf("A=\n");
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
        if(j==0)
            printf("[%f,",X[i][j]);
        else if(j == dim-1)
            printf("%f]\n",X[i][j]);
        else
            printf("%f,",X[i][j]);
        }
    }
}


int diag_domin(float** A, int n)
{
    int i, j;
    float sum;

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

void jacobi(float** A, float* b ,float *x, int n)
{
    int i, j, k, aide;
    float* x_new;
    float diff, sum;

    x_new = (float*)malloc(n * sizeof(float));


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

void gauss_seidel(double** A, double* b, double* x, int dim)
{
    int i, j, k;
    double norme, somme;
    double y[dim];
    for (k = 0; k < 1000; k++)   // 1000 nombre maximal d'iterations
    {
        norme = 0.0;
        for (i = 0; i < dim; i++)
        {
            somme = b[i];
            for (j = 0; j < dim; j++)
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
            printf("Gauss-Seidel a converge en %d iterations\n", k + 1);
            return;
        }
    }

    printf("Gauss-Seidel n'a pas converge en %d iterations\n", k);
}