/*
NAME:		interpolation.c 
DESCRIPTION: 	N-dimensional linear interpolation
AUTHOR:	 	Will Grey
VERSION:	2015-05-05	
LICENSE:	This is free and unencumbered software 
                released into the public domain.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerical.h"

double mult(int dims, double fx, double *pos, double *b){

 int i;
 double val=1.0; 

 for (i = 0; i < dims; i++){
  val *= fabs(1.0-pos[i]-b[i]); 
 } 
 val *= fx;
 
 printf("mult val= %f fx=%f\n",val, fx);
 return val;

}

double interp(int j, int dims, double *pos, double *fx, double *b, int initialise) {
 
 static double fxInterpolated = 0.0;
 static int iter = 0; 
 int i;

 
 if (initialise){
  fxInterpolated = 0.0;
  iter = 0;
 }

 if (j == dims){
  fxInterpolated+=mult(dims,fx[iter],pos,b);
  iter++;
  return fxInterpolated; 
 }

  for (i = 0; i < 2; i++){
   pos[j]=(double)i; 
   interp(j+1,dims,pos,fx,b,0);
  }
 return fxInterpolated; 
}

double linearInterpolation(double *fx, double **a, double *x, int dims){

  double *b, *pos;
  int i;
  double fxInterpolated = 0.0;
 
  b=allocateDoubleVector(dims);
  pos=allocateDoubleVector(dims);
  
  for (i=0; i<dims;i++)
   b[i] = (x[i] - a[i][0]) / (a[i][1] -  a[i][0]);
   
  fxInterpolated=interp(0,dims,pos,fx,b,1); 
  
  return fxInterpolated;
  
}  

double interp2D(int dims, double *fx, double *b){
 
 int i,j;
 double val=0.0, *pos; 
 pos=allocateDoubleVector(dims);
 
  for (i=0; i<2;i++){
   for (j=0; j<2;j++){
    pos[0]=i;
    pos[1]=j;
    val+=mult(dims,fx[j+i*2],pos,b);
   }
  } 
 
 return val;
}

double linearInterpolationQuadratic(double *fx, double x1, double x2, double y1, double y2, double x, double y){

  double b1, b2;
  
  b1= (x - x1)/(x2 - x1);
  b2= (y - y1)/(y2 - y1);
  return (1-b1)*(1-b2)*fx[0] +
         (  b1)*(1-b2)*fx[1] +
         (1-b1)*(  b2)*fx[2] +
         (  b1)*(  b2)*fx[3];

}

double linearInterpolation1D(double f1, double f2,  double x1, double x2, double x){

  double b;
  
  b= (x - x1)/(x2 - x1);
  return (1 - b)*f1 + b*f2;
  
}



