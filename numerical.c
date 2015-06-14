/*
NAME:		numerical.c 
DESCRIPTION: 	For testing the collection of numerical
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
#include <sys/time.h>
#include "numerical.h"

void testUnivariateStats(){
 
 int i,n;
 float v[]={4,11,8,10,6,8};
 unsigned char a[]={4,11,8,10,6,8};
 n=ARRAY_LENGTH(a);
 univariateStats stats;
 
 /*univStats(v, n, &stats);*/

 univStatsByte(a, n, &stats);
 
 printf("\nUnivariate Stats\n");
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

}

void testBivariateStats(){

 
 float x[]={95.0,85.0,80.0,70.0,60.0};
 float y[]={85.0,95.0,70.0,65.0,70.0};
 regressionCoefficients coeffs;
 
 coeffs=regression(x,y,5);
 
 printf("\nBivariate Stats\n"); 
 printf("correlation =%f\n",correlation(x,y,5));
 printf("covariance =%f\n",covariance(x,y,5));
 printf("rmse =%f\n",rmse(x,y,5));
 printf("bias =%f\n",bias(x,y,5));
 printf("m =%f\n",coeffs.m);
 printf("c =%f\n",coeffs.c);
 
}

void testUtility(){

 int i,n;
 short int b[]={4,11,8,10,6,8};
 n=6;
 
 printf("\nUtilities\n");
 printf("Byte swap\n");
 
 byteSwap(b, n, 2);

 for (i=0;i<n;i++){
  printf("%d ",b[i]);
 } 
 printf("\n");

 byteSwap(b, n, 2);
 
 for (i=0;i<n;i++){
  printf("%d ",b[i]);
 } 
 printf("\n");

}

void testSequence(){

 int i;
 long *fib;
 unsigned char *primeNumbmer;
 
 pdf *dist;
 dist=normalDistribution(10.0, 5.0, 4, 0.1, 400);
 
 printf("\nSequence\n");
 printf("Normal distribution\n");
 for (i=0;i<400;i++){
   printf("%f\t%f\n",(dist+i)->x,(dist+i)->fx);
 }  
 
 printf("gamma 1 =%f\n",tgamma(4));
 printf("gamma 1 =%f\n",gammaFunction(4));
 printf("p(a|b) =%f\n",bayes(0.8, 0.01, 0.096));
 printf("factorial =%d\n",factorial(5));
 printf("double factorial =%d\n",doubleFactorial(5));
 
 fib=fibonacci(20);
 primeNumbmer=prime(20);
 
 printf("fibonacci\n");
 for (i=0;i<20;i++){
 printf("%d ",i);
 } 
 printf("\n");

 for (i=0;i<20;i++){
 printf("%ld ",fib[i]);
 } 
 printf("\n");

 printf("prime\n");
 for (i=0;i<20;i++){
  if (primeNumbmer[i] ==1)
   printf("%d ",i);
 } 
  printf("\n");

}



void testVectorAlgebra(){

 double v[]={4,11,8,10};
 double vLength=0;

 double v1[]={3,2,1,-2};
 double v2[]={2,-1,4,1};

 printf("\nMatrix\n");
 vLength=vectorLength(v, 4);
 printf("vector length: %f\n",vLength);

 printf("vector addition\n");
 vectorAddition(v1, v2, 4);
 printVector(v1,4);

}



void testGaussianElimination(){
 
 double **A;
 int n;
 double *x; 
 n=3;
 double b[n];

 A=allocateDoubleMatrix(3, 3);

 A[0][0]=1.0;  A[0][1]=1.0;  A[0][2]=-1.0;
 A[1][0]=1.0;  A[1][1]=2.0;  A[1][2]=1.0;
 A[2][0]=2.0;  A[2][1]=-1.0;  A[2][2]=1.0;
 
 b[0]=2.0;
 b[1]=6.0;
 b[2]=1.0;
 
 x=gaussianElimination (n, A, b); 
 
 fprintf(stdout,"Gaussian Elimination\n");
 printVector(x,n);

}

void testMatrixDeterminant(){

 int n;
 double **C;
 n=3;
 
 C=allocateDoubleMatrix(3, 3);

 C[0][0]=6.0;  C[0][1]=1.0;  C[0][2]=1.0;
 C[1][0]=4.0;  C[1][1]=-2.0;  C[1][2]=5.0;
 C[2][0]=2.0;  C[2][1]=8.0; C[2][2]=7.0;
 
 fprintf(stdout,"Matrix Determinant\n");
 fprintf(stdout,"%f \n", matrixDeterminant(C,n));

}

void testMatrixInverse(){
 
 double **A, **B;
 int n;
 n=3;
 
 A=allocateDoubleMatrix(3, 3);

 A[0][0]=1.0;  A[0][1]=2.0;  A[0][2]=3.0;
 A[1][0]=0.0;  A[1][1]=4.0;  A[1][2]=5.0;
 A[2][0]=1.0;  A[2][1]=0.0;  A[2][2]=6.0;
 
 
 B=matrixInverse(A,n);
 
 fprintf(stdout,"Matrix Inversion\n");
 printMatrix(B,n,n);
 
}

void testLUFactorisation(){
 
 double **A, **L, **U;
 int n;
 n=3;
 
 A=allocateDoubleMatrix(n, n);
 L=allocateDoubleMatrix(n, n);
 U=allocateDoubleMatrix(n, n);

 A[0][0]=6.0;  A[0][1]=0.0;  A[0][2]=2.0;
 A[1][0]=24.0;  A[1][1]=1.0;  A[1][2]=8.0;
 A[2][0]=-12.0;  A[2][1]=1.0;  A[2][2]=-3.0;
 
 LUFactotisation(n,A,U,L);
 
 fprintf(stdout,"Matrix Factorisation U:\n");
 printMatrix(U,n,n);

 fprintf(stdout,"Matrix Factorisation L:\n");
 printMatrix(L,n,n);
 
}

void testLUSolver(){
 
 double **A, *x;
 int n;
 n=3;
 double b[n];
 
 A=allocateDoubleMatrix(n, n);
 

 A[0][0]=6.0;  A[0][1]=0.0;  A[0][2]=2.0;
 A[1][0]=24.0;  A[1][1]=1.0;  A[1][2]=8.0;
 A[2][0]=-12.0;  A[2][1]=1.0;  A[2][2]=-3.0;
 
 b[0]=4.0;
 b[1]=19.0;
 b[2]=-6.0;

 x=LUSolver(A,b,n);
 
 fprintf(stdout,"LU solver U:\n");
 printVector(x,n);

}

void  testInterpolation(){

 double *fx, **a, *x;
 int dims=2;

 x=allocateDoubleVector(dims);
 a=allocateDoubleMatrix(dims,2);
 fx=allocateDoubleVector(dims*2);
 
 x[0]=0.25;
 x[1]=0.4;

 a[0][0]=0;
 a[0][1]=1;
 a[1][0]=0;
 a[1][1]=1;

 fx[0]=1;
 fx[1]=3;
 fx[2]=2;
 fx[3]=4;
 
 printf("\nInterpolation\n");
 printf("%f\n",linearInterpolation(fx, a, x, dims));

} 

void testIntegration(){
 
 long int r;
 r=1867;
 
 printf("\nIntergration\n");
 Integration(funcTest,simpsonsRule,0.0,1.0);
 Integration(funcTest,trapizoidalRule,0.0,1.0);
 printf("monteCarloIntegration1D: %f\n",monteCarloIntegration1D(funcTest,0.0,1.0,1e7));
 printf("simpsonsRule: %f\n",simpsonsRule(funcTest,0.0,1.0,100));
 printf("random01:%f\n",random01());
 printf("random01Self: %f\n",random01Self(&r));
 printf("trapizoidalRule2D: %f\n",trapizoidalRule2D(funcTest5, 0, 1, 0, 1, 10,10)); 
 printf("trapizoidalRule2D: %f\n",trapizoidalRule2D(funcTest6, 1, 2, 3, 5, 40,40));
 
} 

void testRoot(){

 printf("\nRoot finding and minimisation\n");
 golden(functionTest3,0,2);
 golden(functionTest1,3,7);
 golden(functionTest2,5,7);
 golden(poly,-0.5,1.0); 

 brute (poly,-1.0,1.0, 0.0001); 
 brute (functionTest1,3,7, 0.0001); 
 secant(poly,0.0,1.0);
 newton(poly, dpoly, 0.0); 
 regulaFalsi(poly,-1.0,1.0); 
 bisect(poly,-1.0,1.0);
 fixedPointIteration(cosine,1.0); 
 squareRoot(20); 
 
} 

void testTimeStep(){
 
 printf("\nTime Step\n");
 printf("rungeKuttaSecondOrder: %f\n",rungeKuttaSecondOrder(funcTimeStep,5,0.1));
 printf("rungeKuttaSecondOrder: %f\n",rungeKuttaFourthOrder(funcTimeStep,5,0.1));
 printf("matsuno: %f\n",matsuno(funcTimeStep,5,0.1));
 printf("leapFrog: %f\n",leapFrog(funcTimeStep,5,0.1));
 printf("euler: %f\n",euler(funcTimeStep,5,0.1));
 
} 

void testRSignificanceTest(){
 
 printf("\nSignificance Test\n");
 correlationProbability(0.21, 23);
 correlationProbability(0.12, 43);
 correlationProbability(0.4, 30);
 correlationProbability(0.7, 50);

}


void testChiSquared(){

 printf("\nChi Squared test\n");
 printf("chisq: the result is significant at p <%4.3f\n",chiSquaredDistribution(16,3));
 
 double **observed, **observed2;

 observed=allocateDoubleMatrix(3,3);

 observed[0][0]=49;  observed[0][1]=50; observed[0][2]=69; 
 observed[1][0]=24;  observed[1][1]=36; observed[1][2]=38; 
 observed[2][0]=19;  observed[2][1]=22; observed[2][2]=28;

 printf("chisq: the result is significant at p <%4.3f\n",chiProbability(observed,3,3));
 
 observed2=allocateDoubleMatrix(2,2);

 observed2[0][0]=25;  observed2[0][1]=6; 
 observed2[1][0]=8;  observed2[1][1]=15; 
 
 printf("chisq: the result is significant at p <%4.3f\n",chiProbability(observed2,2,2));

}

void testMultipleRegression(){

 int i, n=3, m=6;
 double **X, *y, *x;
 
 X=allocateDoubleMatrix(m, n); 
 y=allocateDoubleVector(m);
 
 X[0][0]=1, X[1][0]=1, X[2][0]=1, X[3][0]=1, X[4][0]=1, X[5][0]=1; 
 X[0][1]=0, X[1][1]=2, X[2][1]=2.5, X[3][1]=1, X[4][1]=4, X[5][1]=7; 
 X[0][2]=0, X[1][2]=1, X[2][2]=2, X[3][2]=3, X[4][2]=6, X[5][2]=2; 
 
 y[0]=5, y[1]=10, y[2]=9, y[3]=0, y[4]=3, y[5]=27; 
 
 x=multipleLinearRegression(X,y,n,m);

  fprintf(stdout,"\nMultiple Regression\n");
 for(i=0;i<n;i++) fprintf(stdout,"%f ",x[i]);
 fprintf(stdout,"\n");
 
 
}

void testSimplex(){
 
 double **x;
 int n=2; 
 
 x = allocateDoubleMatrix(n+1,n);
 x[0][0]=0.0; x[0][1]=0.0;
 x[1][0]=1.2; x[1][1]=0.0;
 x[2][0]=0.0; x[2][1]=0.8;
 
 printf("\nSimplex\n");
 simplexDownhill(testFunctionSimplex,x,n);
 printf("%f %f\n",x[0][0],x[0][1]);

}


int main(int argc, char *argv[])
{

 testSequence();
 testUnivariateStats();
 testBivariateStats();
 testUtility();
 testVectorAlgebra(); 
 testMatrixDeterminant();
 testMatrixInverse();
 testLUFactorisation();
 testLUSolver();
 testGaussianElimination();
 testInterpolation();
 testIntegration();
 testRoot();
 testTimeStep();
 testRSignificanceTest();
 testChiSquared();
 testMultipleRegression();
 testSimplex();

 return (EXIT_SUCCESS);

} 
 



