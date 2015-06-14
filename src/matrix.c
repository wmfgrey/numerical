/*
NAME:		matrix.c 
DESCRIPTION: 	Performs numerous matrix calculations
AUTHOR:	 	Will Grey
VERSION:	2015-05-05	
LICENSE:	This is free and unencumbered software 
                released into the public domain.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>
#include "numerical.h"

double multipleCorrelation(double **X, double *y, int n, int m, double *b){

 int i,j; 
 double *yprime;
 yprime=allocateDoubleVector(m);

 for (i=0; i<m; i++){
  yprime[i]=0.0; 
  for (j=0; j<n; j++){
   yprime[i]+=b[j]*X[i][j];
  }
 }
  return correlation ((float*)y,(float*)yprime,m);
}

double * multipleLinearRegression(double **X, double *y, int n, int m){

 double **A, **XT, **B, *b, *x, **Y;
 int i,j; 

 Y=vectorToMatrix(y,m,1);

 x=allocateDoubleVector(n);
 XT=trasposeMatrix(X,n,m); 
/* XT.X.x = XT.b */
 A=multiplyMatrix(XT,X,m,n,n);
 B=multiplyMatrix(XT,Y,m,n,1); 
 b=matrixToVector(B,n,1);

 for (i=0; i<n; i++){ 
  for (j=0; j<n; j++){
   fprintf(stdout,"%f ",A[i][j]); 
  }
  fprintf(stdout,"\n"); 
 }  
 
 for(i=0;i<n;i++) fprintf(stdout,"%f ",b[i]);
  fprintf(stdout,"\n"); 
 
 x = gaussianElimination (n, A, b); 

 return x;

}

double * LUSolver(double **A, double *b, int n){
 
 double **L, **U, *c, *x;
 
 L=allocateDoubleMatrix(n, n);
 U=allocateDoubleMatrix(n, n);

 LUFactotisation(n,A,U,L);
 c=forwardSubstitution(n, L, b);
 x=backSubstitution(n, U, c);
 
 return x;
}

int LUFactotisation (int n, double **A, double **U, double **L){
 
 int i, j, k;
 double m;
 
 resetMartixDoubleZero(U,n,n);
 resetMartixDoubleZero(L,n,n);

 for (i=0;i<n;i++){
  U[0][i]=A[0][i];
  L[i][i]=1.0;
 }

  for (k = 0; k<n-1; k++){
  for (i = k+1; i<n; i++){ 
  
   m = A[i][k] / A[k][k];
   L[i][k]=m;
   
   for (j=k+1; j<n; j++) 
    U[i][j] = A[i][j] = A[i][j] - m * A[k][j];
   
  }
 }

 for (k = 0; k<n; k++)
  for (i = 0; i<n; i++)
   if(k<i) U[i][k]=0.0; 
 
 return (EXIT_SUCCESS);

}


double ** matrixInverse(double **A, int n){

 double **B, **C;
 double invDet, det;
 
 det=matrixDeterminant(A,n);
 
 if (det!=0){
  invDet=1.0/det;
  B=cofactorMatrix(A,n);
  C=trasposeMatrix(B, n, n);
  B=scalarMultiplyMatrix(C,invDet,n,n); 
 }
 else{
  B=allocateDoubleMatrix(n, n);
  resetMartixDoubleZero(B,n,n);
 }

 return B;
}

double ** cofactorMatrix(double **A, int n){

 double **B;
 int i, j;
 B=allocateDoubleMatrix(n, n);
 
 for (j=0; j<n; j++){
  for (i=0; i<n; i++){
   B[j][i]=pow(-1,i+j)*matrixDeterminant(minorMatrix(A,n,i,j),n-1);
   
  }
 }
 return B;
}


double ** minorMatrix(double **A, int n, int p, int q){

 double **B;
 int k, j, m=0, l=0;
 
 B=allocateDoubleMatrix(n-1, n-1);

 for(j=0;j<n;j++){
  m=0;
  for(k=0;k<n;k++){
     if(k!=p && j!=q){
      B[l][m]=A[j][k];
      m++;
     }
   }
   if(j!=q) l++; 
  }

 return B;
}

double matrixDeterminant(double **A, int n){

 int i;
 double det=0;
 double **B;

 B=allocateDoubleMatrix(n-1, n-1);
 
 if(n==2){
  det=(A[0][0]*A[1][1])-(A[0][1]*A[1][0]);
 }

 else{
  for(i=0;i<n;i++){        
    B=minorMatrix(A, n, i, 0);
    det += pow(-1,i)*A[0][i]*matrixDeterminant(B,n-1);   
   }
  }

 return(det);
}




double * gaussianElimination (int n, double **A, double *b){
 
 int i, j, k;
 double m, *x; 

 for (k = 0; k<n-1; k++){
  for (i = k+1; i<n; i++){ 
  
   m = A[i][k] / A[k][k];
   
   for (j=k+1; j<n; j++) 
    A[i][j] = A[i][j] - m * A[k][j];
   
   b[i] = b[i] - m * b[k];   
  }
 }

 x=backSubstitution (n, A, b);
 
 return x;

}

double * backSubstitution (int n, double **A, double *b){
 
 int i, j;
 double *x;
  
 x=allocateDoubleVector(n);

 x[n-1] = b[n-1]/A[n-1][n-1];

 for(i=n-2; i>=0; i--){

  x[i] = b[i];
  for (j=i+1; j<n; j++) 
   x[i] = x[i] - A[i][j] * x[j];
  
  x[i] = x[i]/A[i][i];

 } 
  
 return x;

}

double * forwardSubstitution (int n, double **A, double *b){
 
 int i, j;
 double *x, a;
  
 x=allocateDoubleVector(n);
  
 for(i=0; i<n; i++){

  a = b[i];
  for (j=0; j<=i; j++) 
   a = a - A[i][j] * x[j];
   
  x[i] = a/A[i][i];
  
 } 
 
 return x;

}





double matrixTrace(double **A, int n){
  
  double sum=0;
  int i;
  for (i=0; i<n; i++)
  sum+=A[i][i];

  return sum;
}

int printMatrix(double **A, int n, int m){
 int i, j;

 for (i=0; i<n; i++){
  for (j=0; j<m; j++){
   fprintf(stdout,"%6.2f\t",A[i][j]);
  }
  fprintf(stdout,"\n");
 }

 return (EXIT_SUCCESS);
}

int printVector(double *a, int n){
 int i;

 for (i=0; i<n; i++){
   fprintf(stdout,"%6.2f\t",a[i]);
 }
 fprintf(stdout,"\n");
 
 return (EXIT_SUCCESS);
}



double ** copyMatrix(double **A, int n, int m){
 
 int i, j;
 double **B;
 B=allocateDoubleMatrix(n, m);

 for (i=0; i<n; i++)
  for (j=0; j<m; j++)
   B[i][j]=A[i][j];
 
  return B;
}

double * copyVector(double *a, int n){
 
 int i;
 double *b;
 b=allocateDoubleVector(n);

 for (i=0; i<n; i++)
   b[i]=a[i];
 
  return b;
}



double * vectorInnerProduct(double *a, double *b, int n)
{

 int i;
 double *c;
 c = allocateDoubleVector(n);
  
 for(i=0;i<n;i++)
  c[i]=a[i]*b[i];
 
 return c;

}

double * vectorDivision(double *a, double *b, int n)
{

 int i;
 double *c;
 c = allocateDoubleVector(n);
  
 for(i=0;i<n;i++)
  c[i]=a[i]/b[i];
 
 return c;

}

double * vectorSubtract(double *a, double *b, int n)
{

 int i;
 double *c;
 c = allocateDoubleVector(n);
  
 for(i=0;i<n;i++)
  c[i]=a[i]-b[i];
 
 return c;

}

double * vectorAddition(double *a, double *b, int n)
{

 int i;
 double *c;
 c = allocateDoubleVector(n);
  
 for(i=0;i<n;i++)
  c[i]=a[i]+b[i];
 
 return c;

}

double * scalarMultiplyVector(double *a, double scalar, int n)
{

 int i;
 double *b;
 b = allocateDoubleVector(n);
 
 for(i=0;i<n;i++)
  b[i]=scalar*a[i];
 
 return b;

}

double * normalVector(double *a, int n)
{

 int i;
 double vLength=0;
 double *b;
 b = allocateDoubleVector(n);

 vLength=vectorLength(a,n);

 for(i=0;i<n;i++)
  vLength=vLength+(a[i]*a[i]);
 
 vLength=sqrt(vLength);

 for(i=0;i<n;i++)
  b[i]=a[i]/vLength;


 return b;

}




double ** addMatrix(double **A, double **B,  int n, int m){
 
 int i, j;
 double **X;
 X=allocateDoubleMatrix(n, m);

 for (i=0; i<n; i++)
  for (j=0; j<m; j++)
   X[i][j]=A[i][j]+B[i][j];
 
  return X;
}

double ** subtractMatrix(double **A, double **B,  int n, int m){
 
 int i, j;
 double **X;
 X=allocateDoubleMatrix(n, m);

 for (i=0; i<n; i++)
  for (j=0; j<m; j++)
   X[i][j]=A[i][j]-B[i][j];
 
  return X;
}

double ** vectorToMatrix(double *b, int n, int m){
 
 int i, j;
 double **A;
 A=allocateDoubleMatrix(n,m);

 for (i=0; i<m; i++)
  for (j=0; j<n; j++)
   A[j][i]=b[j+i*m];

  return A;
}


double *matrixToVector(double **A, int n, int m){
 
 int i, j;
 double *b;
 b=allocateDoubleVector(n*m);

 for (i=0; i<m; i++)
  for (j=0; j<n; j++)
   b[j+i*m]=A[j][i];

  return b;
}

double ** multiplyMatrix(double **A, double **B,  int m, int n1, int n2){
 
 int i, j, k;
 double **X;
 X=allocateDoubleMatrix(n1, n2);

 for (i=0; i<n1; i++)
  for (j=0; j<n2; j++)
   X[i][j]=0.0;

 
 for (k=0; k<n2; k++)
  for (i=0; i<n1; i++)
   for (j=0; j<m; j++)
    X[i][k]+=A[i][j]*B[j][k];
  
  return X;
}


double ** trasposeMatrix(double **A, int n, int m){
 
 int i, j;
 double **AT;
 AT=allocateDoubleMatrix(n, m);

 for (i=0; i<n; i++)
  for (j=0; j<m; j++)
   AT[i][j]=A[j][i];

  return AT;
}


double ** diagonalMatrix(double **A, int n){

 int i, j;
 double **B;
 B=allocateDoubleMatrix(n, n);

 for (i=0;i<n;i++){
  for (j=0;j<n;j++){
   if (i !=j){
    B[i][j]=0.0;
   }
   else{
    B[i][j]=A[i][j];
   }
  }
 }
 return B;
}

double ** scalarMultiplyMatrix(double **A, double scalar, int n, int m){

 int i, j;
 double **B;
 B=allocateDoubleMatrix(n, m);

 for (i=0;i<n;i++)
  for (j=0;j<m;j++)
   B[i][j]=scalar*A[i][j];

 return B;
}


double ** createIdentityMatrix(int n){

 int i, j;
 double **A;

 A=allocateDoubleMatrix(n,n);

 for (i=0;i<n;i++)
  for (j=0;j<n;j++)
   A[i][j]=0.0;

 for (i=0;i<n;i++)
  A[i][i]=1.0;

 return A;
}

int resetMartixDoubleZero(double **A, int n, int m){

 int i, j;

 for (i=0;i<n;i++)
  for (j=0;j<m;j++)
   A[i][j]=0.0;

 return (EXIT_SUCCESS);
}



int resetVectorDoubleZero(double *a, int n){

 int i;

 for (i=0;i<n;i++)
   a[i]=0.0;

 return (EXIT_SUCCESS);
}


double vectorLength(double *a, int n)
{

 int i;
 double vLength=0;
 for(i=0;i<n;i++)
  vLength=vLength+(a[i]*a[i]);
 
 return sqrt(vLength);

}

double * allocateDoubleVector(int n){
 
 int i;
 double *vector;
 if((vector = (double *) calloc(n,sizeof(double)))==NULL) 
  memoryCheck();

 for (i=0;i<n;i++)
   vector[i]=0.0;
 
 return vector;
}


float * allocateFloatVector(int n){

 int i;
 float *vector;
 if((vector = (float *) calloc(n,sizeof(float)))==NULL) 
  memoryCheck();
 
 for (i=0;i<n;i++)
   vector[i]=0.0;

 return vector;
}


short int * allocateShortVector(int n){
 
 int i;
 short int *vector;
 if((vector = (short int *) calloc(n,sizeof(short int)))==NULL) 
  memoryCheck(); 

 for (i=0;i<n;i++)
   vector[i]=0.0;
 
 return vector;
}


long int * allocateLongVector(int n){

 int i;
 long int *vector;
 if((vector = (long int *) calloc(n,sizeof(long int)))==NULL) 
  memoryCheck();

 for (i=0;i<n;i++)
   vector[i]=0.0;
 
 return vector;
}


unsigned char * allocateCharVector(int n){

 int i;
 unsigned char *vector;
 if((vector = (unsigned char *) calloc(n,sizeof(unsigned char)))==NULL) 
  memoryCheck(); 

  for (i=0;i<n;i++)
   vector[i]=0.0;
 
 return vector;
}

void freeDoubleVector(double * vector){ free(vector);}
void freeFloatVector(float * vector){ free(vector); }
void freeCharVector(unsigned char * vector){ free(vector); }
void freeLongVector(long int * vector){ free(vector);}
void freeShortVector(short int * vector){ free(vector);}

double ** allocateDoubleMatrix(int i, int j){
 
 int k,n,m;
 double ** matrix;

 if((matrix = (double **) calloc(i, sizeof(double *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateDoubleVector(j);
 
 for (n=0;n<i;n++)
  for (m=0;m<j;m++)
   matrix[n][m]=0;

 return matrix;
}

float ** allocateFloatMatrix(int i, int j){
 
 int k,n,m;
 float ** matrix;

 if((matrix = (float **) calloc(i, sizeof(float *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateFloatVector(j);

 for (n=0;n<i;n++)
  for (m=0;m<j;m++)
   matrix[n][m]=0;

 
 return matrix;
}

unsigned char ** allocateCharMatrix(int i, int j){
 
 int k,n,m;
 unsigned char ** matrix;

 if((matrix = (unsigned char **) calloc(i, sizeof(unsigned char *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateCharVector(j);

 for (n=0;n<i;n++)
  for (m=0;m<j;m++)
   matrix[n][m]=0;
 
 return matrix;
}

long int ** allocateLongMatrix(int i, int j){
 
 int k,n,m;
 long int ** matrix;

 if((matrix = (long int **) calloc(i, sizeof(long int *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateLongVector(j);

 for (n=0;n<i;n++)
  for (m=0;m<j;m++)
   matrix[n][m]=0;
 
 return matrix;
}

short int ** allocateShortMatrix(int i, int j){
 
 int k,n,m;
 short int ** matrix;

 if((matrix = (short int **) calloc(i, sizeof(short int *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateShortVector(j);

 for (n=0;n<i;n++)
  for (m=0;m<j;m++)
   matrix[n][m]=0;
 
 return matrix;
}


void freeDoubleMatrix(double ** matrix, int rows){
 int i;
 for (i =0; i < rows; i++) freeDoubleVector(matrix[i]);
 free(matrix);
}


void freeFloatMatrix(float ** matrix, int rows){
 int i;
 for (i =0; i < rows; i++) freeFloatVector(matrix[i]);
 free(matrix);
}


void freeLongMatrix(long int ** matrix, int rows){
 int i;
 for (i =0; i < rows; i++) freeLongVector(matrix[i]);
 free(matrix);
}

void freeShortMatrix(short int ** matrix, int rows){
 int i;
 for (i =0; i < rows; i++) freeShortVector(matrix[i]);
 free(matrix);
}

void freeCharMatrix(unsigned char ** matrix, int rows){
 int i;
 for (i =0; i < rows; i++) freeCharVector(matrix[i]);
 free(matrix);
}

