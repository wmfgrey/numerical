/*
NAME:		simplex.c 
DESCRIPTION: 	Using Nelder Mead for N-dimensional
                function mininisation.
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
#include "numerical.h"

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

