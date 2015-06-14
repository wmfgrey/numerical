/* 
NAME:		bivariateStats.c 
DESCRIPTION: 	Calculating bivariate statistics including regression,
                correlation, covariance, rmse and bias.  
AUTHOR:	 	Will Grey
VERSION:	2015-05-15
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

regressionCoefficients regression(float *xData, float *yData, int n){
	
 regressionCoefficients coeffs;
 int i;
 double SumX=0,SumY=0,SumX2=0,SumXY=0;

 for(i=0; i<n; i++){
															
  SumX += (double)xData[i];
  SumY += (double)yData[i];
  SumX2 += ((double)xData[i] * (double)xData[i]);
  SumXY += ((double)xData[i] * (double)yData[i]); 
 
 }
 


 coeffs.m = ((n*SumXY)-(SumY*SumX)) / ((n*SumX2)-(SumX*SumX));
 coeffs.c = (SumY - (coeffs.m * SumX)) / n;

 return coeffs;

}

double correlation (float *xData, float *yData, int n){

 int i;
 double sumProduct=0, xSumSquares=0, ySumSquares=0;
 double xMean, yMean;

 xMean = mean(xData,n);
 yMean = mean(yData,n);

 for(i=0; i<n; i++){
  
  sumProduct += ((double)xData[i] - xMean) * ((double)yData[i] - yMean);
  xSumSquares += pow((double)xData[i] - xMean,2);
  ySumSquares += pow((double)yData[i] - yMean,2);

 }
 printf("%f %f %f %f %f\n",sumProduct, xSumSquares, ySumSquares, xMean, yMean);
 
 return  SQUARE(sumProduct / ( sqrt(xSumSquares * ySumSquares) ) );

}

double covariance (float *xData, float *yData, int n){

 int i;
 double sumProduct=0;
 double xMean, yMean;

 xMean = mean(xData,n);
 yMean = mean(yData,n);

 for(i=0; i<n; i++)
  
  sumProduct += ((double)xData[i] - xMean) * ((double)yData[i] - yMean);
  
 return sumProduct / (n-1);

}

double rmse (float *xData, float *yData, int n){

 int i;
 double sumSquaredError=0;
 
 for(i=0; i<n; i++)
  
  sumSquaredError += SQUARE((double)xData[i] -(double)yData[i]);
  
 return sqrt(sumSquaredError/n);

}

double bias (float *xData, float *yData, int n){

 double xMean, yMean;

 xMean = mean(xData,n);
 yMean = mean(yData,n);
  
 return xMean-yMean;

}




