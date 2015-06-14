/* 
NAME:	univariateStats.c 
DESCRIPTION: 	Calculating univariate statistics.  
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

int univStats(float *data, int n, univariateStats *stats){
	 
 stats->mean=mean(data,n);
 stats->median=median(data,n);
 stats->max=max(data,n); 
 stats->min=min(data,n); 
 stats->stdDev=stdDev(data,n);
 stats->variance=variance(data,n); 
 stats->skewness=skewness(data,n); 
 stats->kurtosis=kurtosis(data,n);
 stats->range=range(data,n);  
 stats->upperQuartile=upperQuartile(data,n); 
 stats->lowerQuartile=lowerQuartile(data,n);
 stats->interQuartileRange=interQuartileRange(data,n); 
 stats->bestEstimateOfStDev=bestEstimateOfStDev(data,n);
 stats->StandardErrorOfStDev=StandardErrorOfStDev(data,n);
 stats->StandardErrorOfMean=StandardErrorOfMean(data,n); 

 return EXIT_SUCCESS;

}

int univStatsByte(unsigned char *data, int n, univariateStats *stats){
 
 float * fData;
 stats->mode=mode(data,n);
 fData=convertByte2Float(data,n);
 univStats(fData,n,stats);

 return EXIT_SUCCESS;

}


unsigned char mode (unsigned char *data, int n)
{

 int i;
 long *hist, maxval=0;
 unsigned char mode=0;
 
 hist=histogram (data,n);

  for(i=0;i<256;i++) {	
   if(hist[i] > maxval){
    maxval = hist[i];
    mode = i;
   }
  } 
		
 return mode;

}


long * histogram (unsigned char *data, int n)
{

 int i;
 long *hist;
 if((hist  = (long *) calloc(256, sizeof(long)))== NULL) memoryCheck();
 
for(i=0;i<256;i++)	
  hist[i] = 0;

  for (i=0;i<n;i++)
    hist[*(data+i)]++;
		
 return hist;

}

double variance (float *data, int n)
{
  return SQUARE((double)stdDev(data,n));
}

double bestEstimateOfStDev(float *data, int n){
 
 return (double)(stdDev(data,n)/SQUARE(((double)n/((double)n-1))));

}

double StandardErrorOfStDev(float *data, int n){
 
 return (double)(bestEstimateOfStDev(data,n)/SQUARE(2*n));

}

double StandardErrorOfMean(float *data, int n){
 
 return (double)(bestEstimateOfStDev(data,n)/sqrt(n));

}

double interQuartileRange (float *data, int n)
{
  return (double)(upperQuartile(data,n)-lowerQuartile(data,n));
}

double lowerQuartile (float *data, int n)
{
  int lq; 	
  sort(data,n);
	
  if((double)((0.25*(n+1)-(int)(0.25*(n+1))>=0.5))){
    lq=ceil(0.25*(n + 1));
  }
  else{
    lq=floor(0.25*(n + 1));
  }
		
  return (double)*(data + lq - 1);

}

double upperQuartile (float *data, int n)
{
  int uq; 	
  sort(data,n);
	
  if((double)((0.75*(n+1)-(int)(0.75*(n+1))>=0.5))){
    uq=ceil(0.75*(n + 1));
  }
  else{
    uq=floor(0.75*(n + 1));
  }
		
  return (double)*(data + uq - 1);

}

double median (float *data, int n)
{
  int val; 
  double median;	
  sort(data,n);
	
  if((n % 2) == 0){
   val = (int)n / 2;
   median = (double)(((double)*(data+val)+(double)*(data+val-1)) / 2);
  }
  
  if((n % 2) != 0){
   val = (int)n / 2;
   median = (double)*(data + val);
  }

  return median;

}

int sort(float *data, int n)
{

 int i, sorted=0, check=1; 
 float holdValue;
	
 while( sorted == 0 ){

  check=1;
		
  for (i=1; i < n; i++){
		
   if (data[i] < data[i - 1])	{
		
     holdValue = data[i];
     data[i] = data[i - 1];
     data[i - 1] = holdValue;	
     check=0;
    }	
  }
		
   if(check == 1) sorted = 1;
 }
	
 return (EXIT_SUCCESS);

}


double range (float *data, int n){
 
 return (double)max(data,n)-(double)min(data,n); 

}

double stdDev (float *data, int n){
 
 double sumSquares=0, avg=0;
 int i;
 
 avg=mean(data, n);

 for (i=0;i<n;i++){
  sumSquares += SQUARE((double)*(data+i)-avg);
 }

 return sqrt(sumSquares/n); 

}

double skewness (float *data, int n){
 
 double sumSquares=0, sumSkewness=0, avg=0;
 int i;
 
 avg=mean(data, n);

 for (i=0;i<n;i++){
  sumSquares += SQUARE((double)*(data+i)-avg);
  sumSkewness += pow((double)*(data+i)-avg,3);
 }
 
 return (double)(sumSkewness/n) / (pow((double)(sumSquares/n),3)); 

}

double kurtosis (float *data, int n){
 
 double sumSquares=0, sumKurtosis=0, avg=0;
 int i;
 
 avg=mean(data, n);

 for (i=0;i<n;i++){
  sumSquares += SQUARE((double)*(data+i)-avg);
  sumKurtosis += pow((double)*(data+i)-avg,4);
 }
 
 return (double)(sumKurtosis/n) / (pow((double)(sumSquares/n),4)); 

}


double mean (float *data, int n){
 
 double sum=0;
 int i;

 for (i=0;i<n;i++){
  sum = sum+(double)*(data+i);
 }

 return sum/n; 

}

double sum (float *data, int n){
 
 double sum=0;
 int i;

 for (i=0;i<n;i++){
  sum = sum+(double)*(data+i);
 }

 return sum; 

}


double max (float *data, int n){
 
 double m=-1e12;
 int i;

 for (i=0;i<n;i++){
  m = MAXIMUM(m,(double)*(data+i));
 }

 return m; 

}

double min (float *data, int n){
 
 double m=1e12;
 int i;

 for (i=0;i<n;i++){
  m = MINIMUM(m,(double)*(data+i));
 }

 return m; 

}


