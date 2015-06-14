#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define ERR_TOL    1e-7
#define MAX_ITERS  200
#define TRIALS     1e5
#define SEED       1e3  

double gaussianIntegration(double (*f)(double), double, double);
double monteCarloIntegration1D(double (*f)(double), double, double, int);
double monteCarloIntegration(double (*f)(double *, int), double *, double *, int, int);
double trapizoidalRule(double (*f)(double), double, double, int);
double simpsonsRule(double (*f)(double), double, double, int);
double Integration(double (*f)(double), double (*method)(double (*f)(double), double, double, int), double, double);
double  trapizoidalRule2D(double (*f)(double, double), double, double, double, double, int, int);
double funcTest(double);
double funcTest5(double, double);
double funcTest6(double, double);
double random01();
double random01Self(long int *r);

int main(int argc, char *argv[])
{
 long int r;
 r=1867;

 Integration(funcTest,simpsonsRule,0.0,1.0);
 Integration(funcTest,trapizoidalRule,0.0,1.0);
 printf("%f\n",monteCarloIntegration1D(funcTest,0.0,1.0,1e7));
 printf("%f\n",simpsonsRule(funcTest,0.0,1.0,100));
 printf("%f\n",gaussianIntegration(funcTest,0.0,1.0));
 printf("%f\n",random01());
 printf("%f\n",random01Self(&r));
 printf("%f\n",trapizoidalRule2D(funcTest5, 0, 1, 0, 1, 10,10)); 
 printf("%f\n",trapizoidalRule2D(funcTest6, 1, 2, 3, 5, 40,40));
 
 return (EXIT_SUCCESS);

} 

double monteCarloIntegration(double (*f)(double *, int), double *a, double *b, int n, int dims){

 long int i, sum=0;
 double volume=1.0, *x;
 
 for (i=0;i<dims;i++)
  volume *= b[i]-a[i];
 
 srand48(SEED); 

 x = (double *) calloc(dims, sizeof(double));

 for (i=0;i<n;i++){
  
  for(i=0;i<dims;i++)
   x[i]=a[i] + drand48() * (b[i]-a[i]);
  
 sum+=f(x,dims);

 }
  
 return (sum / (double)n) * volume;

}

double monteCarloIntegration1D(double (*f)(double), double a, double b, int n){

 long int i;
 double x, sum=0;

 srand48(SEED); 

 for (i=0;i<n;i++){

  x=drand48();
  sum+=f(x);
  
 }
  
 return (sum / (double)n) * (b-a);

}



double  simpsonsRule(double (*f)(double), double a, double b, int n)
{  
 
 int i;
 double sum=0, h,integral, x;
 
 sum = (f(a) + f(b))/3.0;

 h=(b-a)/(double)n;
 
 for (i=1;i<n;i++){
  x = a+(i*h);
  if (pow(-1.0,i) < 0){
   sum+=f(x)*(4.0/3.0);
  }
  else{
   sum+=f(x)*(2.0/3.0);
  }  
 }

 integral = h*sum;
 return integral;
   
} 



double  trapizoidalRule(double (*f)(double), double a, double b, int n)
{  
 
 int i;
 double sum=0, h,integral, x;
 
 sum=(f(a) + f(b))/2;
 h=(b-a)/(double)n;
 
 for (i=1;i<n;i++){
  x = a+(i*h);
  sum+=f(x);
 }

 integral = h*sum;
 return integral;
   
}

double  trapizoidalRule2D(double (*f)(double, double), double a, double b, double c, double d, int n, int m)
{  
 
 int i, j;
 double sum[6], h, k, integral=0, x, y;
  
 h=(b-a)/(double)m;
 k=(d-c)/(double)n;
 
 for (i=0;i<6;i++)
  sum[i]=0;

 sum[0]=f(a,c) + f(b,c) + f(a,d) + f(b,d);
 
 for (i=1;i<m;i++){
  x = a+(i*h);
  sum[1]+=f(x,c);
  sum[2]+=f(x,d);
 }

 for (j=1;j<n;j++){
  y = c+(j*k);
  sum[3]+=f(a,y);
  sum[4]+=f(b,y);
 }
 
 for (i=1;i<5;i++)
  sum[i] *=2;
 
 for (i=1;i<m;i++){
  for (j=1;j<n;j++){
    x = a+(i*h);
    y = c+(j*k);
    sum[5]+=f(x,y);
  }
 }
 
 sum[5] *=4;
 
 for (i=0;i<6;i++)
  integral +=sum[i];

 integral = ((h*k)/4 ) *  integral;
 
 return integral;
   
} 


 

double Integration(double (*f)(double), double (*method)(double (*f)(double), double, double, int), double a, double b)
{  
 
 int i=2;
 double integral1=0.0, integral2=0.0;
 
 integral1=method(f,a,b,1);

 while ((fabs(integral1-integral2) > ERR_TOL) && (i < MAX_ITERS)){
  i++;
  integral1=integral2;
  integral2=method(f,a,b,i);
  
 }
 fprintf(stdout,"Approximation to integral %f\n",integral2); 
 if (i>=MAX_ITERS) fprintf(stdout,"Desired accuracy not achieved\n"); 

 return integral2;
   
} 

double gaussianIntegration(double (*f)(double), double a, double b){
 
 double x1, x2, I=0;
 x1 = (b+a-(b-a))/sqrt(3)/2.0;
 x2 = (b+a+(b-a))/sqrt(3)/2.0;
 I=(b-a)*(f(x1)+f(x2))/2;
 return I;

}

double funcTest(double x){
 
 return sqrt((x*x)+1); 

}

double random01(){
 
 return (double) ((double)rand() / (double)RAND_MAX+1.0);

}

double random01Self(long int *r){
 
 int a=106,c=1283,m=6075;
 *r = (*r*a+c)%m;
 return (double)*r/(double)m;
 
}

double funcTest5(double x, double y){
 
 return 8 * exp(-pow(x,2)-pow(y,4)); 

}

double funcTest6(double x1, double x2){
 
 return x1 * pow(x2,3); 

}



