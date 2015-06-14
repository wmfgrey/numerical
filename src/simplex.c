
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define ERR_TOL    1e-6
#define MAX_ITERS  200


int memoryCheck(void);	
void freeDoubleVector(double * );
void freeFloatVector(float * );
void freeCharVector(unsigned char * );
void freeLongVector(long int * );
void freeShortVector(short int * );
double ** allocateDoubleMatrix(int, int);
float ** allocateFloatMatrix(int, int);
unsigned char ** allocateCharMatrix(int, int);
short int ** allocateShortMatrix(int, int);
long int ** allocateLongMatrix(int, int);
void freeDoubleMatrix(double **, int);
void freeFloatMatrix(float **, int);
void freeLongMatrix(long int **, int);
void freeShortMatrix(short int **, int);
void freeCharMatrix(unsigned char **, int);
double * allocateDoubleVector(int);
float * allocateFloatVector(int);
short int * allocateShortVector(int);
long int * allocateLongVector(int);
unsigned char * allocateCharVector(int);

int repositionSimplex(short int *, double **, int);
int sortSimplex(double *, int, short int *);
int calcCentroid(double **, int, double *);
double testFunctionSimplex(double *, int);
double simplexDownhill (double (*f)(double *, int), double **, int);


void testSimplex(){
 
 double **x;
 int n=2; 
 
 x = allocateDoubleMatrix(n+1,n);
 x[0][0]=0.0; x[0][1]=0.0;
 x[1][0]=1.2; x[1][1]=0.0;
 x[2][0]=0.0; x[2][1]=0.8;

 simplexDownhill(testFunctionSimplex,x,n);
 printf("%f %f\n",x[0][0],x[0][1]);

}

int main(int argc, char *argv[])
{

 testSimplex();

 return (EXIT_SUCCESS);

} 


double simplexDownhill (double (*f)(double *, int), double **x, int n){
 
 double *xm, *xr, *xe, *xoc, *xic, *fx, tol=1e6, fr, fe, foc, fic;
 double alpha=1, beta=2, gamma=-0.5, delta=0.5;
 int i, j, iters=0;
 short int *xpos;

 xpos = allocateShortVector(n+1);
 fx = allocateDoubleVector(n+1);
 xm = allocateDoubleVector(n);
 xr = allocateDoubleVector(n);
 xe = allocateDoubleVector(n);
 xoc = allocateDoubleVector(n);
 xic = allocateDoubleVector(n);
 
 while (iters < MAX_ITERS && tol > ERR_TOL){

 for (i=0;i<n+1;i++){
  fx[i] = f(x[i],n);
  xpos[i]=i; 
 }
 
 sortSimplex(fx,n+1,xpos);
 repositionSimplex(xpos, x, n);
 calcCentroid(x,n,xm);
 
 for (i=0;i<n+1;i++){
  for (j=0;j<n;j++){  
  printf("%f ",x[i][j]); 
  }
  printf("%f\t",fx[i]); 
 }
 printf("\n"); 
 
 /* reflection */ 

 for (i=0;i<n;i++)
  xr[i] = xm[i] + alpha * (xm[i] - x[n][i]);
 
 fr=f(xr,n);
 
 if (fx[0] <= fr && fr < fx[n-1]){
   for (i=0;i<n;i++) 
    x[n][i]=xr[i];  
 }

 else if (fr < fx[0]){
  
  /* expansion */
  
  for (i=0;i<n;i++)
   xe[i] = xm[i] + beta * (xr[i] - xm[i]);
  fe=f(xe,n); 
  
  if (fe < fr){
    for (i=0;i<n;i++) 
     x[n][i]=xe[i];
   }

   else{
    for (i=0;i<n;i++) 
     x[n][i]=xr[i];
   } 

 }

 
 else if (fx[n-1] <= fr && fr <  fx[n]){
  
  /* outside contraction */
  
  for (i=0;i<n;i++)
   xoc[i] = xm[i] + gamma * (x[n][i] - xm[i]);
  
  foc =f(xoc,n);

  if (foc <= fr){
   for (i=0;i<n;i++) 
     x[n][i]=xoc[i];
  }

  else{
   /* reduce */ 
   for (j=1;j<n+1;j++) 
    for (i=0;i<n;i++) 
     x[j][i] = x[0][i] + delta * (x[j][i]-x[0][i]);

  }
 

 }

 else if (fr >= fx[n]){

  /* inside contraction */

  for (i=0;i<n;i++) 
   xic[i] = xm[i] - gamma * (x[n][i] - xm[i]);
  
  fic =f(xic,n);

  if (fic < fx[n]){
   for (i=0;i<n;i++) 
     x[n][i]=xic[i];
  } 
  
  else{
   /* reduce */ 
   for (j=1;j<n+1;j++) 
    for (i=0;i<n;i++) 
     x[j][i] = x[0][i] + delta * (x[j][i]-x[0][i]);
  }
 }

  tol = fabs(fx[0]-fx[n]);
  iters++;

}


 freeShortVector(xpos);  
 freeDoubleVector(xe); 
 freeDoubleVector(xoc); 
 freeDoubleVector(xic);

 return fx[0];
}


int calcCentroid(double **x, int n, double *xm){

 int i,j;
 for (i=0;i<n;i++) xm[i]=0;
 
 for (i=0;i<n;i++){
  for (j=0;j<n;j++){
   xm[j]+=x[i][j];
  }
  xm[j]=xm[j]/n;
 }

 return (EXIT_SUCCESS);

}


int repositionSimplex(short int *xpos, double **x, int n){

 int i,j;
 double **xorder;
 xorder = allocateDoubleMatrix(n+1,n);

 for (i=0;i<n+1;i++){
  for (j=0;j<n;j++){  
   xorder[i][j]=x[xpos[i]][j];
  }
 }
 
 for (i=0;i<n+1;i++)
  for (j=0;j<n;j++)
   x[i][j]=xorder[i][j];
 

 freeDoubleMatrix(xorder, n);


 return (EXIT_SUCCESS);
}



int sortSimplex(double *data, int n, short int *xpos)
{

 int i, sorted=0, check=1; 
 double holdValue;
 	
 while( sorted == 0 ){

  check=1;
		
  for (i=1; i < n; i++){
		
   if (data[i] < data[i - 1])	{
		
     holdValue = data[i];
     data[i] = data[i - 1];
     data[i - 1] = holdValue;	
     
     holdValue = xpos[i];
     xpos[i] = xpos[i - 1];
     xpos[i - 1] = holdValue;	
     
     check=0;
    }	
  }
		
   if(check == 1) sorted = 1;
 }
	
 return (EXIT_SUCCESS);

}

double testFunctionSimplex(double *x,int n){

 return pow(x[0],2) - (4 * x[0]) + pow(x[1],2) - x[1] - (x[0] * x[1]);
  
}



/* utility routines */




double * allocateDoubleVector(int i){
 double *vector;
 if((vector = (double *) calloc(i,sizeof(double)))==NULL) 
  memoryCheck();  
 return vector;
}


float * allocateFloatVector(int i){
 float *vector;
 if((vector = (float *) calloc(i,sizeof(float)))==NULL) 
  memoryCheck();  
 return vector;
}


short int * allocateShortVector(int i){
 short int *vector;
 if((vector = (short int *) calloc(i,sizeof(short int)))==NULL) 
  memoryCheck();  
 return vector;
}


long int * allocateLongVector(int i){
 long int *vector;
 if((vector = (long int *) calloc(i,sizeof(long int)))==NULL) 
  memoryCheck();  
 return vector;
}


unsigned char * allocateCharVector(int i){
 unsigned char *vector;
 if((vector = (unsigned char *) calloc(i,sizeof(unsigned char)))==NULL) 
  memoryCheck();  
 return vector;
}

void freeDoubleVector(double * vector){ free(vector);}
void freeFloatVector(float * vector){ free(vector); }
void freeCharVector(unsigned char * vector){ free(vector); }
void freeLongVector(long int * vector){ free(vector);}
void freeShortVector(short int * vector){ free(vector);}

double ** allocateDoubleMatrix(int i, int j){
 
 int k;
 double ** matrix;

 if((matrix = (double **) calloc(i, sizeof(double *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateDoubleVector(j);
 
 return matrix;
}

float ** allocateFloatMatrix(int i, int j){
 
 int k;
 float ** matrix;

 if((matrix = (float **) calloc(i, sizeof(float *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateFloatVector(j);
 
 return matrix;
}

unsigned char ** allocateCharMatrix(int i, int j){
 
 int k;
 unsigned char ** matrix;

 if((matrix = (unsigned char **) calloc(i, sizeof(unsigned char *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateCharVector(j);
 
 return matrix;
}

long int ** allocateLongMatrix(int i, int j){
 
 int k;
 long int ** matrix;

 if((matrix = (long int **) calloc(i, sizeof(long int *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateLongVector(j);
 
 return matrix;
}

short int ** allocateShortMatrix(int i, int j){
 
 int k;
 short int ** matrix;

 if((matrix = (short int **) calloc(i, sizeof(short int *)))==NULL)
  memoryCheck();

 for (k=0; k< i; k++)matrix[k] = allocateShortVector(j);
 
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

int memoryCheck(void)	
{	
 fprintf(stderr,"\nUnable to allocate memory.\n\n");
 exit(EXIT_FAILURE);
}


