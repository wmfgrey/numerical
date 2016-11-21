#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "numerical.h"


double golden (double (*f)(double), double a, double b){
 
 double goldenRatio = (sqrt(5)-1)/2;
 int i=0;
 double x1, x2;
 
 x1 = b - goldenRatio * (b - a);
 x2 = a + goldenRatio * (b - a);

 while((fabs(b - a) > ERR_TOL) && (i < MAX_ITERS)){
  
  if (f(x2) > f(x1)){
    b=x2;
    x2=x1;
    x1 = b - goldenRatio * (b - a);
  }
  else{
    a=x1;
    x1=x2;
    x2 = a + goldenRatio * (b - a); 
  }
   i++;

  }

 fprintf(stdout,"Golden: %f %f %d\n", f(fabs(b+a)/2), fabs(b+a)/2, i); 
 return fabs(b+a)/2;

}

double goldenFit (double a, double b, float *x, float *y, int num, double (*f)(float, float)){
 
 double goldenRatio = (sqrt(5)-1)/2;
 int i=0;
 double x1, x2;
 
 x1 = b - goldenRatio * (b - a);
 x2 = a + goldenRatio * (b - a);

 while((fabs(b - a) > ERR_TOL) && (i < MAX_ITERS)){
  
  if (merit_function_golden(&x2,x,y,num,f) > merit_function_golden(&x1,x,y,num,f)){
    b=x2;     
    x2=x1;
    x1 = b - goldenRatio * (b - a);
  }
  else{
    a=x1;
    x1=x2;
    x2 = a + goldenRatio * (b - a); 
  }
   i++;

  }
 
 x2=fabs(b+a)/2;

 fprintf(stdout,"Golden: %f %f %d\n", merit_function_golden(&x2,x,y,num,f), fabs(b+a)/2, i); 
 return fabs(b+a)/2;

}

float merit_function_golden(double *a, float *x, float *y, int n, double (*f)(float, float)){

 float e=0.0;
 int i;
 
 for (i=0;i<n;i++)
   e+=pow(y[i]-f(*a,x[i]),2);  

 return e;

}




