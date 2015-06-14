/*
NAME:		timeStep.c 
DESCRIPTION: 	Collection of time stepping integration routines.	
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "numerical.h"


double euler(double (*f)(double, double), int timeSteps, double h){
 
 double y=1.0, x;
 int i;

 for (i=0;i<timeSteps;i++){
  x=h*i;
  y=  y + h * f(x,y);

 }
 
 return y;

}

double rungeKuttaSecondOrder(double (*f)(double, double), int timeSteps, double h){
 
 double y=1.0, x, k1, k2;
 int i;

 for (i=0;i<timeSteps;i++){
  
  x = h * i;
  k1 = h * f(x,y);
  k2 = h * f(x+h,y+k1);
  y = y + ((k1 + k2)/2);
  
 }
 
 return y;

}

double leapFrog(double (*f)(double, double), int timeSteps, double h){
 
 double y, y1, x=h, y0=1.0;
 int i;

 y1=y0+h*f(0.0,y0);

 for (i=1;i<timeSteps;i++){
  y =  y0 + h * 2 * f(x,y1);
  x=h*i;
  y0 = y1;
  y1=y;
 }
 
 return y;

}
 
double matsuno(double (*f)(double, double), int timeSteps, double h){
 
 double y=1.0, x, yprime;
 int i;

 for (i=0;i<timeSteps;i++){
  x=h*i;
  yprime = y + h * f(x,y);
  y = y + f(x,yprime) * h;

 }
 
 return y;

}

double rungeKuttaFourthOrder(double (*f)(double, double), int timeSteps, double h){

 double k1, k2, k3, k4, y=1.0, x;
 int i;
 
 for (i=0;i<timeSteps;i++){
  x=h*i;
  k1 = f(x,y) * h;
  x=x+h/2;  
  k2 = f(x,y+k1/2) * h;
  k3 = f(x,y+k2/2) * h;
  x=x+h/2;
  k4 = f(x,y+k3) * h;
  y = y + ((k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0);
}
 return y;
}









