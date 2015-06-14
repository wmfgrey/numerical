/*Multiple Regression */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double multipleCorrelation(double **, double *, int, int, double *);
double * multipleLinearRegression(double **, double *, int, int);
double * gaussianElimination (int, double **, double *);
double * backSubstitution (int, double **, double *);
double ** allocateDoubleMatrix(int, int);
double * allocateDoubleVector(int);
double ** multiplyMatrix(double **, double **, int, int, int);
double ** trasposeMatrix(double **, int, int);
double * matrixToVector(double **, int, int);
double ** vectorToMatrix(double *, int, int);
double correlation (float *, float *, int);
double mean (float *, int);

int main(int argc, char *argv[]){ 

 int i, n=3, m=6;
 double **X, *y, *x;
 
 X=allocateDoubleMatrix(m, n); 
 y=allocateDoubleVector(m);
 
 X[0][0]=1, X[1][0]=1, X[2][0]=1, X[3][0]=1, X[4][0]=1, X[5][0]=1; 
 X[0][1]=0, X[1][1]=2, X[2][1]=2.5, X[3][1]=1, X[4][1]=4, X[5][1]=7; 
 X[0][2]=0, X[1][2]=1, X[2][2]=2, X[3][2]=3, X[4][2]=6, X[5][2]=2; 
 
 y[0]=5, y[1]=10, y[2]=9, y[3]=0, y[4]=3, y[5]=27; 
 
 x=multipleLinearRegression(X,y,n,m);

 for(i=0;i<n;i++) fprintf(stdout,"%f ",x[i]);
 fprintf(stdout,"\n");
 
 return EXIT_SUCCESS;
}

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

double correlation (float *xData, float *yData, int n){

 int i;
 double sumProduct=0, xSumSquares=0, ySumSquares=0;
 double xMean, yMean;

 xMean = mean(xData,n);
 yMean = mean(yData,n);

 for(i=0; i<n; i++){
  
  sumProduct += ((double)xData[i] - xMean) * ((double)yData[i] - yMean);
  xSumSquares += pow((double)xData[i] - xMean,2);
  ySumSquares += pow((double)yData[i] - yMean,2);

 }
 printf("%f %f %f %f %f\n",sumProduct, xSumSquares, ySumSquares, xMean, yMean);
 
 return  pow(sumProduct / ( sqrt(xSumSquares * ySumSquares) ),2);

}

double mean (float *data, int n){
 
 double sum=0;
 int i;

 for (i=0;i<n;i++){
  sum = sum+(double)*(data+i);
 }

 return sum/n; 

}



