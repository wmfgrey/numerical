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

double test_linear_model(float a, float h){
 
 return a * h;

}

double test_linear_model_2(double *p, double h){

 return p[0] +  p[1] * h;

}

double testFunctionSimplexFit(double *x){

double e;

 e = pow((x[0] + (1.707106 * x[1]) - 12.500000),2) + 
     pow((x[0] + (2.490711 * x[1]) - 4.166667),2) +
     pow((x[0] + (3.605550 * x[1]) - 31.250000),2) +
     pow((x[0] + (4.297620 * x[1]) - 81.250000),2) +
     pow((x[0] + (5.656853 * x[1]) - 112.500000),2);

 return e;
}



double testFunctionSimplexRosenbrock(double *x){

 return 100 * pow((x[1] - pow(x[0],2)),2) + pow(1-x[0],2);
  
}


double testFunctionSimplex(double *x){

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




