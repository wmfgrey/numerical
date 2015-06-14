/*
NAME:		testFunctions.c 
DESCRIPTION: 	Range of mathematical functions for testing
                the numerical routines.	
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

double funcTimeStep(double x, double y){
 
 return x + y; 

}

double testFunctionSimplex(double *x,int n){

 return pow(x[0],2) - (4 * x[0]) + pow(x[1],2) - x[1] - (x[0] * x[1]);
  
}

/* Root */

double functionTest1(double x){

 return -1.0*pow(x-4.1,2);

}

double functionTest2(double x){

 return x * sin(x) - 2 * cos(x);

}

double functionTest3(double x){

 return x*x + x - 2*sqrt(x);

}

double poly(double x)
{
 
 double a=3.0, b=1.0, c=-5.0, d=-1.0;
 return d + x*(c + x*(b + x*a));

}

double dpoly(double x)
{
 
 double a=3.0, b=1.0, c=-5.0;
 return c+ x*(2 * b + x * 3 * a);

}
  
double cosine(double x){

 x=cos(x);
 return x;

}

double exponent(double x){

 x=exp(-1.0*x);
 return x;

}

/* Integration */

double funcTest(double x){
 
 return sqrt((x*x)+1); 

}

double funcTest5(double x, double y){
 
 return 8 * exp(-pow(x,2)-pow(y,4)); 

}

double funcTest6(double x1, double x2){
 
 return x1 * pow(x2,3); 

}




