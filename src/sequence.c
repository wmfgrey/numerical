/*
NAME:		sequence.c 
DESCRIPTION: 	Sequencing and other routines.	
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

#define BIG 2e6      
#define SEED 1e3  

int combination(int n, int r){
 return factorial(n)/(factorial(r)*factorial(n-r));;
}

int factorial(int n)
{	
 int f=0;
 
 if( n >= 0)  f=1;
 
 while (n > 1){
   f *= n;
   n-=1;   
  }
 	
  return f;
}

int doubleFactorial(int n)
{	
 int f=0;
 
 if( n >= 0 || n == -1)  f=1;
 
 while (n > 1){
   f *= n;
   n-=2;   
  }
 	
  return f;
}
      
double gammaFunction(double k){
/*Closed form of Gamma function because 
  we are only dealing with integers */
 k=k*2;

 
 return (double)doubleFactorial(k-2) * sqrt(PI) / pow(2,(k-1)/2);

 


}

double incompleteGammaFunction(double chi, int n){
 
/* http://mathworld.wolfram.com/Chi-SquaredDistribution.html */

 double sum=0;
 int k;
 
 for (k=0;k<n;k++)
  sum += pow(chi,k) / factorial(k);

 return factorial(n-1) * (1 - exp(-1.0*chi) * sum);

}

unsigned char *prime( int n){

 int i,j, remainder;
 unsigned char *primeNumber;
 if((primeNumber  = (unsigned char *) calloc(n, sizeof(unsigned char)))== NULL) memoryCheck();
 
 primeNumber[0]=0;
 primeNumber[1]=1;
 
 for (j=2;j<=n;j++){
   
  primeNumber[j]=1;
 
  for (i=2;i<j;i++){
   remainder=j % i;
   
   if ((remainder == 0)){
    primeNumber[j]=0;
    break;   
   }
  }
 }

 return primeNumber;

}


long *fibonacci(int maxFibonacciNumber){
 
 int i, num1=0, num2=1; 
 long *fibonacciNumber;
 if((fibonacciNumber  = ( long*) calloc(maxFibonacciNumber, sizeof(long)))== NULL) memoryCheck();
 
 fibonacciNumber[0]=num1;
 fibonacciNumber[1]=num2; 
 
 for (i=2;i<=maxFibonacciNumber;i++){
  fibonacciNumber[i] = num1 + num2;
  num1=num2;
  num2=fibonacciNumber[i];

 }
  
 return fibonacciNumber;
}

double bayes(double probBA, double probA, double probBnotA){

 double probNotA, probB, probAB; 

 probNotA = 1.0-probA;

 probB = ( probBA * probA ) + ( probBnotA * probNotA);
 probAB = ( probBA * probA ) / probB;

 return probAB;
}


pdf *normalDistribution(double mean, double stdDev, int numStdDev, double binSize, int numBins){

 pdf *dist;
 int i; 
 double x;
 
 if((dist = (pdf *) calloc(numBins, sizeof(pdf)))== NULL) memoryCheck();
 
  for (i=0;i<numBins;i++){
  x = mean - (stdDev * (double)numStdDev) + ((double)i * binSize); 
  (dist+i)->x = x;
  (dist+i)->fx = ((1/sqrt(2 * PI * SQUARE(stdDev)) * exp(-1.0*(SQUARE(x-mean)/(2*SQUARE(stdDev))))) * binSize);
 }
 return dist;
}












