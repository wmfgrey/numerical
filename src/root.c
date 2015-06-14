/*
NAME:		root.c 
DESCRIPTION: 	collection of tools for find roots and
                minimum of functions.
                routines.	
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

double secant(double (*f)(double), double x1, double x2)
{  
 
 /* Secant method */
 
 double  f1, f2=1.0, x;
 int i=0;
 
 while((fabs(f2) > ERR_TOL) && (i < MAX_ITERS)){
  f1=f(x1);
  f2=f(x2);
  x=x2-((f2*(x2-x1))/(f2-f1));
  x1=x2;
  x2=x;
  i++;    
 }   
 
 fprintf(stdout,"Secant: %f %f %d\n",f1, x, i); 
 return x;
   
} 

double newton(double (*f)(double), double (*fd)(double), double x)
{  
 
 /* Newton Raphson method */
 
 double fx, dfx;
 int i;
 
 i=0;
 fx=1.0;

 while((fabs(fx) > ERR_TOL) && (i < MAX_ITERS)){
  fx=f(x);
  dfx=fd(x);
  x=x-fx/dfx;
  i++; 
    
 }   
 
 fprintf(stdout,"Newton: %f %f %d\n",fx, x, i); 
 return x;
   
} 


double regulaFalsi(double (*f)(double), double a, double b)
{  
 
 /* Regula falsi */
 
 double  fa, fb, x, fx;
 int i;
 
 i=0;
 fa=f(a);
 fb=f(b);
 fx=1.0;

 while((fabs(fx) > ERR_TOL) && (i < MAX_ITERS)){
  
  x=((a*fb) - (b*fa)) / (fb-fa); 
  fx=f(x);
  if ((fa*fx) > 0.0){
   a=x;
   fa=fx;
  }
  else{
   b=x;
   fb=fx;
  } 
  i++; 
    
 }   
 
 fprintf(stdout,"Regula Falsi: %f %f %d\n",fx, x, i);  
 return x;
   
} 


double bisect(double (*f)(double), double a, double b)
{  
 
 /* bisection method */
 
 double fa, x, fx;
 int i;
 
 i=0;
 fa=f(a);
 fx=1.0;

 while((fabs(fx) > ERR_TOL) && (i < MAX_ITERS)){
  x=0.5*(a+b);
  fx=f(x);
  if ((fa*fx) > 0.0){
   a=x;
   fa=fx;
  }
  else{
   b=x;
  } 
  i++;       
 }   
 
 fprintf(stdout,"Bisection: %f %f %d\n",fx, x, i);
 return x;
   
} 


double brute (double (*f)(double), double a, double b, double eps){

 int i, n;
 double x, xmin=1e6, fx, fxmin=1e6;
 
 n=(int)((b-a)/eps);
 
 for (i=0;i<n;i++){
  x=a+(i*eps);
  fx=f(x);
 
  if (fabs(fx) < fxmin){
   fxmin=fabs(fx);
   xmin=x;
  }
   
 }
 
 fprintf(stdout,"Brute: %f %f %d\n",fxmin, xmin, i);
 return xmin;

}


double fixedPointIteration(double (*f)(double), double x)
{  
 
 /* fixed point iteration 
  Rewrite f(x)=0 in the form of g(x)=x
  Iterate x(m+1)=g(x(m)) */
 
 int i;
 double  err;
 
 i=0;
 err=1.0;

 while((err > ERR_TOL) && (i < MAX_ITERS)){
  x=f(x);
  err=fabs(f(x)-x);  
  i++;
} 
  
 fprintf(stdout,"Fixed Point Iteration: %f %f %d\n",x, err, i);   
 return x;
   
} 


double squareRoot(double a){

 return newtonSqaureRoot(sqr, dsqr, 1.0, a);
  
}

double sqr(double x, double a){
 return x*x-a;

}

double dsqr(double x){
 
 return 2*x;

}

double newtonSqaureRoot(double (*f)(double, double), double (*fd)(double), double x, double a)
{  
 
 /* Newton Raphson method */
  
 double fx, dfx;
 int i;
 
 i=0;
 fx=1.0;

 while((fabs(fx) > ERR_TOL) && (i < MAX_ITERS)){
  fx=f(x,a);
  dfx=fd(x);
  x=x-fx/dfx;
  i++; 
    
 }   
 
 fprintf(stdout,"Newton sqrt: %f %f %d\n",x, fx, i); 
 return x;
   
} 


 

