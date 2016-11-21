/*
NAME:		simplex.c 
DESCRIPTION: 	Using Nelder Mead for N-dimensional
                function mininisation.
AUTHOR:	 	Will Grey
xeRSION:	2015-05-05	
LICENSE:	This is free and unencumbered software 
                released into the public domain.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "numerical.h"


double simplex_fit (double **x, int n, double *k, double *y, int num, double (*f)(double *, double)){
 
 double *xm, *xr, *xe, *xc, *fx, tol=1e6, fr, fe, fc, fmin, m;
 double alpha=1, beta=2, gamma=0.5, delta=0.5;
 int i, j, iters=0;
 short int *xpos;
  
 xpos = allocateShortVector(n+1);
 fx = allocateDoubleVector(n+1);
 xm = allocateDoubleVector(n);
 xr = allocateDoubleVector(n);
 xe = allocateDoubleVector(n);
 xc = allocateDoubleVector(n);
 
 while (iters < MAX_ITERS && tol > ERR_TOL_FIT){

 for (i=0;i<n+1;i++){
  fx[i] = merit_function(x[i],k,y,num,f);
  xpos[i]=i; 
 }
 
 sortSimplex(fx,n+1,xpos);
 repositionSimplex(xpos, x, n);
 calcCentroid(x,n,xm);
 
 fprintf(stdout,"Iteration %d: \n",iters); 
 for (i=0;i<n+1;i++){
  for (j=0;j<n;j++){  
  fprintf(stdout,"%f ",x[i][j]); 
  }
  fprintf(stdout,"%f %f\n",fx[i], tol); 
 }
 

 /* reflection */ 

 for (i=0;i<n;i++)
  xr[i] = xm[i] + alpha * (xm[i] - x[n][i]); 
 fr=merit_function(xr,k,y,num,f);
  
 if (fx[0] <= fr && fr < fx[n-1]){
   for (i=0;i<n;i++) 
    x[n][i]=xr[i];  
 }

 else if (fr < fx[0]){
  
  /* expansion */
  
  for (i=0;i<n;i++)
   xe[i] = xm[i] + beta * (xr[i] - xm[i]);
  fe=merit_function(xe,k,y,num,f); 
  
  if (fe < fx[n-1]){
    for (i=0;i<n;i++) 
     x[n][i]=xe[i];
   }

   else{
    for (i=0;i<n;i++) 
     x[n][i]=xr[i];
   } 
 }


 else{
  
  /* contraction */
  
  for (i=0;i<n;i++)
   xc[i] = xm[i] + gamma * (x[n][i] - xm[i]);
  
  fc=merit_function(xc,k,y,num,f);

  if (fc <= fr){
   for (i=0;i<n;i++) 
     x[n][i]=xc[i];
  }

  else{
   /* reduce */ 
   for (j=1;j<n+1;j++) 
    for (i=0;i<n;i++) 
     x[j][i] = x[0][i] + delta * (x[j][i]-x[0][i]);

  } 
 }

  m = 0.0;
  for (i=0;i<n+1;i++)
   m += fx[i];
  
  m = m/(n+1);

  tol=0.0;
  for (i=0;i<n+1;i++)
   tol +=pow((fx[i]-m),2.0)/(n);

  tol = sqrt(tol);
  iters++;

 }

 fmin=fx[0];

 fprintf(stdout,"Minimum\n"); 
 for (i=0;i<n;i++) 
  fprintf(stdout,"%f ",x[0][i]);
 fprintf(stdout,"%f\n",fmin); 
  
 fprintf(stdout,"Iterations: %d\n",iters);

 freeShortVector(xpos);  
 freeDoubleVector(xe); 
 freeDoubleVector(xc); 
 freeDoubleVector(xr);
 freeDoubleVector(xm);
 freeDoubleVector(fx);

 return fmin;
}

double merit_function(double *x, double *k, double *y, int n, double (*f)(double *, double)){

 double e=0.0;
 int i;
 
 for (i=0;i<n;i++){
  e+=pow(y[i]-f(x,k[i]),2);
 }  
 return e;

}

double simplex (double (*f)(double *), double **x, int n){
 
 double *xm, *xr, *xe, *xc, *fx, tol=1e9, fr, fe, fc, fmin, m;
 double alpha=1, beta=2, gamma=0.5, delta=0.5;
 int i, j, iters=0;
 short int *xpos;
  
 xpos = allocateShortVector(n+1);
 fx = allocateDoubleVector(n+1);
 xm = allocateDoubleVector(n);
 xr = allocateDoubleVector(n);
 xe = allocateDoubleVector(n);
 xc = allocateDoubleVector(n);
 
 while (iters < MAX_ITERS && tol > ERR_TOL){

 for (i=0;i<n+1;i++){
  fx[i] = f(x[i]);
  xpos[i]=i; 
 }
 
 sortSimplex(fx,n+1,xpos);
 repositionSimplex(xpos, x, n);
 calcCentroid(x,n,xm);
 
 fprintf(stdout,"Iteration %d: \n",iters); 
 for (i=0;i<n+1;i++){
  for (j=0;j<n;j++){  
  fprintf(stdout,"%f ",x[i][j]); 
  }
  fprintf(stdout,"%f %f\n",fx[i], tol); 
 }
 

 /* reflection */ 

 for (i=0;i<n;i++)
  xr[i] = xm[i] + alpha * (xm[i] - x[n][i]); 
 fr=f(xr);
  
 if (fx[0] <= fr && fr < fx[n-1]){
   for (i=0;i<n;i++) 
    x[n][i]=xr[i];  
 }

 else if (fr < fx[0]){
  
  /* expansion */
  
  for (i=0;i<n;i++)
   xe[i] = xm[i] + beta * (xr[i] - xm[i]);
  fe=f(xe); 
  
  if (fe < fx[n-1]){
    for (i=0;i<n;i++) 
     x[n][i]=xe[i];
   }

   else{
    for (i=0;i<n;i++) 
     x[n][i]=xr[i];
   } 
 }


 else{
  
  /* contraction */
  
  for (i=0;i<n;i++)
   xc[i] = xm[i] + gamma * (x[n][i] - xm[i]);
  
  fc =f(xc);

  if (fc <= fr){
   for (i=0;i<n;i++) 
     x[n][i]=xc[i];
  }

  else{
   /* reduce */ 
   for (j=1;j<n+1;j++) 
    for (i=0;i<n;i++) 
     x[j][i] = x[0][i] + delta * (x[j][i]-x[0][i]);

  } 
 }

 m = 0.0;
 for (i=0;i<n+1;i++)
  m += fx[i];
  
 m = m/(n+1);

 tol=0.0;
 for (i=0;i<n+1;i++)
  tol +=pow((fx[i]-m),2.0)/(n);

 tol = sqrt(tol);
 iters++;

}

 fmin=fx[0];

 fprintf(stdout,"Minimum\n"); 
 for (i=0;i<n;i++) 
  fprintf(stdout,"%f ",x[0][i]);
 fprintf(stdout,"%f\n",fmin); 
  
 fprintf(stdout,"Iterations: %d\n",iters);

 freeShortVector(xpos);  
 freeDoubleVector(xe); 
 freeDoubleVector(xc); 
 freeDoubleVector(xr);
 freeDoubleVector(xm);
 freeDoubleVector(fx);

 return fmin;
}

int calcCentroid(double **x, int n, double *xm){

 int i,j;
 
 for (j=0;j<n;j++) {
  xm[j]=0.0;
  for (i=0;i<n;i++) 
   xm[j] += x[i][j];
  xm[j] = xm[j]/n;
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

