/*
 
NAME:		numerical.c 
DESCRIPTION: 	This is code that has a number of numerical routines.	
AUTHOR:	 	Will Grey
VERSION:	2015-05-05
TO COMPILE:	gcc numerical.c -o numerical -Wall -lm
EXAMPLE:  	
LICENSE:	This is free and unencumbered software released into the public domain.
			 				
*/

#include "numerical.h"

int main(int argc, char *argv[])
{
 

 int i,n,j;
 float v[]={4,11,8,10,6,8};
 unsigned char a[]={4,11,8,10,6,8};
 short int b[]={4,11,8,10,6,8};
 univariateStats stats;
 long *fib;
 unsigned char *primeNumbmer;
 float x[]={95.0,85.0,80.0,70.0,60.0};
 float y[]={85.0,95.0,70.0,65.0,70.0};
 regressionCoefficients coeffs;
 
 n=ARRAY_LENGTH(v);
 
 printf("chisq: %4.3f\n ",chiSquaredDistribution(7.63,3));
 
 for (j=0;j<=8;j++){
 for (i=1;i<=10;i++){
 
 printf("%4.3f ",chiSquaredDistribution((double)j,i));
 }
 printf("\n");
}

 printf("gamma 1 =%f\n",tgamma(4));
 printf("gamma 1 =%f\n",gammaFunction(4));

 coeffs=regression(x,y,5);

 printf("correlation =%f\n",correlation(x,y,5));
 printf("covariance =%f\n",covariance(x,y,5));
 printf("rmse =%f\n",rmse(x,y,5));
 printf("bias =%f\n",bias(x,y,5));
 printf("m =%f\n",coeffs.m);
 printf("c =%f\n",coeffs.c);
 

 pdf *dist;
 dist=normalDistribution(10.0, 5.0, 4, 0.1, 400);

 for (i=0;i<400;i++){
   printf("%f\t%f\n",(dist+i)->x,(dist+i)->fx);
 }  
 
 
 printf("p(a|b) =%f\n",bayes(0.8, 0.01, 0.096));
 printf("pi =%f\n",monteCarloIntegrationPi());
 printf("factorial =%d\n",factorial(5));
 printf("double factorial =%d\n",doubleFactorial(5));
 
 fib=fibonacci(20);
 primeNumbmer=prime(20);
 
 printf("fibonacci\n");
 for (i=0;i<20;i++){
 printf("%d %ld \n",i,fib[i]);
 } 

 printf("prime\n");
 for (i=0;i<20;i++){
  if (primeNumbmer[i] ==1)
   printf("%d \n",i);
 } 

 printf("Byte swap\n");
 
 byteSwap(b, n, 2);

 for (i=0;i<n;i++){
  printf("%d\n",b[i]);
 } 

 byteSwap(b, n, 2);
 
 for (i=0;i<n;i++){
  printf("%d\n",b[i]);
 } 

 /*univStats(v, n, &stats);*/

 univStatsByte(a, n, &stats);

 printf("n:\t%d\n",n); 
 printf("Mean:\t%f\n",stats.mean);
 printf("Median:\t%f\n",stats.median);
 printf("Mode:\t%d\n",stats.mode);
 
 printf("Max:\t%f\n",stats.max); 
 printf("Min:\t%f\n",stats.min); 

 printf("Stdev:\t%f\n",stats.stdDev);
 printf("Variance:\t%f\n",stats.variance); 
 printf("Skew:\t%f\n",stats.skewness); 
 printf("Kurt:\t%f\n",stats.kurtosis);
 printf("Range:\t%f\n",stats.range); 
 
 printf("UpperQuartile:\t%f\n",stats.upperQuartile); 
 printf("LowerQuartile:\t%f\n",stats.lowerQuartile);
 printf("InterQuartieRange:\t%f\n",stats.interQuartileRange); 

 printf("BestEstimateOfStDev:\t%f\n",stats.bestEstimateOfStDev);
 printf("StandardErrorOfStDev:\t%f\n",stats.StandardErrorOfStDev);
 printf("StandardErrorOfMean:\t%f\n",stats.StandardErrorOfMean);

 printf("sort:\n"); 
 sort(v,n);
 for (i=0;i<n;i++){
  printf("%f\n",v[i]); 
 }

 printf("Mean2:\t%f\n",mean(convertByte2Float(a,n),n));


 return (EXIT_SUCCESS);

} 
 



