#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

int complexConvert(float *, float *, float *, float *, int);
int RealDFT(float *, float *, int);
int ForwardDFT(float *, float *, float *, float *,  int);
int InverseDFT(float *, float *, float *, float *, int);

int main(int argc, char *argv[])
{

 
  return EXIT_SUCCESS;
}

int complexConvert(float *real, float *imag, float *mag, float *phase, int n){
 
 int i;

 for (i=0;i<n;i++){
  *(mag + i)=sqrt(*(real + i) * *(real + i) + *(imag + i) * *(imag + i));
  *(phase + i)=atanf(*(imag + i) / *(real + i));
 }

 return EXIT_SUCCESS;

}

int RealDFT(float *realIn, float *realOut, int n){

  int i, k;
  
 for (i=0;i<n;i++){
  realOut[i]=0.0;
  for (k=0;k<n;k++){
   realOut[k]+=realIn[k]*cos(2*PI*k*i/n);
  }
 }

  return EXIT_SUCCESS;
}

int ForwardDFT(float *realIn, float *imagIn, float *realOut, float *imagOut,  int n){

  int i, k;
  
 for (i=0;i<n;i++){
  realOut[i]=0.0;
  for (k=0;k<n;k++){
   realOut[i]+=realIn[k]*cos(2*PI*k*i/n)+imagIn[k]*sin(2*PI*k*i/n);
   imagOut[i]+=-1.0*realIn[k]*sin(2*PI*k*i/n)+imagIn[k]*cos(2*PI*k*i/n);
  }
 }

  return EXIT_SUCCESS;
}

int InverseDFT(float *realIn, float *imagIn, float *realOut, float *imagOut, int n){

 int i;
  
 ForwardDFT(realIn,imagIn,realOut,imagOut,n);
 
 for (i=0;i<n;i++){
  realOut[i]=realOut[i]/n;
  imagOut[i]=imagOut[i]/n;
 }
 
 return EXIT_SUCCESS;
}



