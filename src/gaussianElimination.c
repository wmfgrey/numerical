#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int gaussianElimination (int, double **, double *, double *);
double ** allocateDoubleMatrix(int, int);
double * allocateDoubleVector(int);

int main(int argc, char *argv[])
{

 int i, n;
 n=3;
 double b[n], x[n];
 
 double **A;

 A=allocateDoubleMatrix(3, 3);
 A[0][0]=1.0;  A[0][1]=1.0;  A[0][2]=-1.0;
 A[1][0]=1.0;  A[1][1]=2.0;  A[1][2]=1.0;
 A[2][0]=2.0;  A[2][1]=-1.0;  A[2][2]=1.0;
 
 b[0]=2.0;
 b[1]=6.0;
 b[2]=1.0;
 
 gaussianElimination (n, A, b, x); 
 
 
 for(i=0;i<n;i++) fprintf(stdout,"%f ",x[i]);
 fprintf(stdout,"\n");

 return EXIT_SUCCESS;
 
}


int gaussianElimination (int n, double **A, double *b, double *x)
{
 
 int i, j, k;
 double m; 
 
 for (k = 0; k<n-1; k++){
  for (i = k+1; i<n; i++){ 
  
   m = A[i][k] / A[k][k];
   
   for (j=k+1; j<n; j++) 
    A[i][j] = A[i][j] - m * A[k][j];
   
   b[i] = b[i] - m * b[k];   
  }
 }
 
 
 x[n-1] = b[n-1]/A[n-1][n-1];

 for(i=n-2; i>=0; i--){

  x[i] = b[i];
  for (j=i+1; j<n; j++) 
   x[i] = x[i] - A[i][j] * x[j];
  
  x[i] = x[i]/A[i][i];

 } 
  
 return EXIT_SUCCESS;

}

double ** allocateDoubleMatrix(int i, int j){
 
 int k;
 double ** matrix;

 matrix = (double **) calloc(i, sizeof(double *));

 for (k=0; k< i; k++)matrix[k] = allocateDoubleVector(j);
 
 return matrix;
}

double * allocateDoubleVector(int i){
 double *vector;
 vector = (double *) calloc(i,sizeof(double));
 return vector;
}

