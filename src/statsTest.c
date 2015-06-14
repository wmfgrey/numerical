#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define T_TABLE_PROBABILITY 4
#define T_TABLE_DOF 36
#define CHI_TABLE_DOF 10
#define CHI_TABLE_PRBABILITY 3


double correlationTTest(double, int);
double rSignifcanceTest(double, int, int);
int correlationProbability(double, int);
double chiSquaredDistribution(double,int);
double chiProbability(double **, int, int);
double chiSquared(double **, int, int);

double ** allocateDoubleMatrix(int, int);
double * allocateDoubleVector(int);

void testRSignificanceTest(){
 
 correlationProbability(0.21, 23);
 correlationProbability(0.12, 43);
 correlationProbability(0.4, 30);
 correlationProbability(0.7, 50);

}


void testChiSquared(){

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

int main(int argc, char *argv[]){

 testRSignificanceTest();
 testChiSquared();

 return (EXIT_SUCCESS);

}

double chiProbability(double **observed, int categories, int dataTypes){
 
 int dof = (categories-1)*(dataTypes-1);
 double chi;
 chi=chiSquared(observed,categories,dataTypes);
 printf("chisq: %4.2f dof: %d\n",chi, dof);
 return chiSquaredDistribution(chi,dof);
 
}

double chiSquared(double **observed, int categories, int dataTypes){
 
 double *observedRow, *observedColumn, **expected, observedTotal=0;
 int i,j; 
 double chi=0;
 
 observedRow=allocateDoubleVector(categories);
 observedColumn=allocateDoubleVector(dataTypes);
 expected=allocateDoubleMatrix(categories,dataTypes);
 
 for (i=0;i<categories;i++){
  for (j=0;j<dataTypes;j++){
   observedRow[i]+=observed[j][i];
   observedColumn[j]+=observed[j][i];
   observedTotal+=observed[j][i];
  } 
 }
 

 for (i=0;i<categories;i++){
  for (j=0;j<dataTypes;j++){
   expected[j][i]=observedRow[i]*observedColumn[j]/observedTotal;
   chi += pow((observed[j][i] - expected[j][i]),2) / expected[j][i];
  }
 } 
  
 return chi;

}


double chiSquaredDistribution(double chi, int dof){
 
 int i; 
 double p=1.0;
 double  chi2Distribution[CHI_TABLE_DOF][CHI_TABLE_PRBABILITY]= 
  {{3.84, 6.64, 10.83},{5.99, 9.21, 13.82},{7.82, 11.34, 16.27},
   {9.49, 13.28, 18.47},{11.07, 15.09, 20.52},{12.59, 16.81, 22.46},
   {14.07, 18.48, 24.32},{15.51, 20.09, 26.12},{16.92, 21.67, 27.88},
   {18.31, 23.21, 29.59}};

 double prob[CHI_TABLE_PRBABILITY]={0.05, 0.01, 0.001};

 if (dof > CHI_TABLE_DOF) dof=CHI_TABLE_DOF;
 
 for (i=0;i<CHI_TABLE_PRBABILITY;i++)
  if (chi >= chi2Distribution[dof-1][i])
   p=prob[i];
   
 return p;
} 


double correlationTTest(double r, int n){

return r * sqrt(n-2) /  sqrt(1-pow(r,2));

}

int correlationProbability(double r, int N){
 
 printf("\nr is: %f\n",r);
 printf("Number of samples is: %d\n",N);
 printf("The t-value is: %f\n",correlationTTest(r,N)); 
 printf("The number of degrees of freedom is: %d\n",N-2);
 printf("For the one sided test the result is significant at p <%f\n",rSignifcanceTest(r,N,0));
 printf("For the two sided test the result is significant at p <%f\n",rSignifcanceTest(r,N,2));

 return (EXIT_SUCCESS);

}



double rSignifcanceTest(double r, int N, int tail){
 
 double prob=1.0,t;
 int i,j,df; 
 
 double tDist[T_TABLE_DOF][T_TABLE_PROBABILITY] =
 {{31.82, 63.66, 318.31, 636.62}, {6.965, 9.925, 22.327, 31.599},
  {4.541, 5.841, 10.215, 12.924}, {3.747, 4.604, 7.173, 8.610},
  {3.365, 4.032, 5.893, 6.869},   {3.143, 3.707, 5.208, 5.959},
  {2.998, 3.499, 4.785, 5.408},   {2.896, 3.355, 4.501, 5.041},
  {2.821, 3.250, 4.297, 4.781},   {2.764, 3.169, 4.144, 4.587},
  {2.718, 3.106, 4.025, 4.437},   {2.681, 3.055, 3.930, 4.318},
  {2.650, 3.012, 3.852, 4.221},   {2.624, 2.977, 3.787, 4.140},
  {2.602, 2.947, 3.733, 4.073},   {2.583, 2.921, 3.686, 4.015},
  {2.567, 2.898, 3.646, 3.965},   {2.552, 2.878, 3.610, 3.922},
  {2.539, 2.861, 3.579, 3.883},   {2.528, 2.845, 3.552, 3.850},
  {2.518, 2.831, 3.527, 3.819},   {2.508, 2.819, 3.505, 3.792},
  {2.500, 2.807, 3.485, 3.768},   {2.492, 2.797, 3.467, 3.745},
  {2.485, 2.787, 3.450, 3.725},   {2.479, 2.779, 3.435, 3.707},
  {2.473, 2.771, 3.421, 3.690},   {2.467, 2.763, 3.408, 3.674},
  {2.462, 2.756, 3.396, 3.659},   {2.457, 2.750, 3.385, 3.646},
  {2.423, 2.704, 3.307, 3.551},   {2.390, 2.660, 3.232, 3.460},
  {2.374, 2.639, 3.195, 3.416},   {2.364, 2.626, 3.174, 3.390},
  {2.330, 2.581, 3.098, 3.300},   {2.326, 2.576, 3.090, 3.291}};

 int dfTable[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                 20,21,22,23,24,25,26,27,28,29,30,40,60,80,100,1000,1e6};

 double p[]={0.01,0.005,0.001,0.0005};
 
 df=N-2;
 
 t=correlationTTest(r,N);

 for (j=1;j<T_TABLE_DOF;j++){
   if ((df >= dfTable[j-1]) && (df < dfTable[j])){
   
    for (i=0;i<T_TABLE_PROBABILITY-1;i++){
     if (t >= tDist[j-1][i] && t < tDist[j-1][i+1])
      prob=p[i];
    } 
    if (t >= tDist[j-1][T_TABLE_PROBABILITY-1]) prob=p[T_TABLE_PROBABILITY-1]; 
   } 
  }

 if (df >= dfTable[T_TABLE_DOF-1]){
  for (i=0;i<T_TABLE_PROBABILITY-1;i++){
     if (t >= tDist[T_TABLE_DOF-1][i] && t < tDist[T_TABLE_DOF-1][i+1])
      prob=p[i];
    } 
    if (t >= tDist[T_TABLE_DOF-1][T_TABLE_PROBABILITY-1]) prob=p[T_TABLE_PROBABILITY-1]; 
  }
 
  if (tail==2) prob*=2.0;
  return prob;    
 }




double ** allocateDoubleMatrix(int i, int j){
 
 int k, l;
 double ** matrix;

 matrix = (double **) calloc(i, sizeof(double *));
 

 for (k=0; k< i; k++)matrix[k] = allocateDoubleVector(j);

 for (k=0;k<i;k++)
  for (l=0;l<j;l++)
   matrix[l][k]=0;

 
 return matrix;
}

double * allocateDoubleVector(int i){
 int j;
 double *vector;
 vector = (double *) calloc(i,sizeof(double)); 
  
 
 for (j=0;j<i;j++)
  vector[j]=0;
   
 return vector;
}


