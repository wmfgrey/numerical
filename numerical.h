/*
NAME:		numerical.h 
DESCRIPTION: 	header file for function prototypes	
AUTHOR:	 	Will Grey
VERSION:	2015-05-05	
LICENSE:	This is free and unencumbered software 
                released into the public domain.
*/


#define MAXIMUM(x,y) ((x) > (y) ? (x) : (y) )
#define MINIMUM(x,y) ((x) < (y) ? (x) : (y) )
#define SQUARE(x) ((x) * (x))
#define ARRAY_LENGTH(x)  (sizeof(x) / sizeof((x)[0]))
#define ABS(x) ((x) >= 0 ? (x) : (-x))

#define PI 3.1415926535897932
#define E 2.7182818284590452
#define MAX_PATH_NAME_LENGTH 1000 
#define NUM_DISTRIBUTION_BINS 100

#define ERR_TOL_FIT 1e-4
#define ERR_TOL    1e-9
#define MAX_ITERS  200
#define TRIALS     1e5
#define SEED       1e3  
#define GOLDEN_RATIO (SQRT(5)-1)/2

#define T_TABLE_PROBABILITY 4
#define T_TABLE_DOF 36
#define CHI_TABLE_DOF 10
#define CHI_TABLE_PRBABILITY 3

/* univariateStats.c */

typedef struct {
 double mean;
 double median;
 unsigned char mode;
 double max;
 double min;
 double range;
 double stdDev;
 double variance;
 double skewness;
 double kurtosis;
 double upperQuartile;
 double lowerQuartile;
 double interQuartileRange;
 double bestEstimateOfStDev;
 double StandardErrorOfStDev;
 double StandardErrorOfMean;
} univariateStats;

typedef struct {
 double fx;
 double x;
} pdf;

typedef struct  {
 double m;
 double c;
} regressionCoefficients;

int sort(float *, int);
double mean(float  *, int);
unsigned char mode (unsigned char *, int);
long * histogram (unsigned char *, int);
double max(float  *, int);
double min(float  *, int);
double range(float  *, int);
double stdDev(float  *, int);
double skewness(float  *, int);
double kurtosis(float  *, int);
double lowerQuartile(float  *, int);
double upperQuartile(float  *, int);
double interQuartileRange(float  *, int);
double median(float  *, int);
double variance (float  *data, int n);
double bestEstimateOfStDev(float  *, int);
double StandardErrorOfStDev(float  *, int);
double StandardErrorOfMean(float  *, int);
int univStats(float  *, int, univariateStats *);
int univStatsByte(unsigned char *, int, univariateStats *);
double sum (float *, int);
int plotHisto(float *, int, float *, float *, float *, int, long *);

/* utility.c */

float * convertShort2Float(unsigned short int *, int);
float * convertLong2Float(unsigned long int *, int);
float * convertDouble2Float(double *, int);
float * convertByte2Float(unsigned char *, int);

int counter(unsigned int, unsigned int);
int createEnviHeader(char *, char *, int, int, int, int);
int byteSwap(void *, int, int);
FILE *openFile(char *, char *, char *);
int memoryCheck(void);
double timer(void);
float *readFloatData(char *, char *, int);

/* sequence.c */

int factorial(int);
int doubleFactorial(int);
unsigned char *prime(int);
long *fibonacci(int);
double bayes(double, double, double);
pdf* normalDistribution(double, double, int, double, int);
double gammaFunction(double k);
int combination(int,int);


/* bivariateStats.c */
regressionCoefficients regression(float *, float *, int);
double correlation (float *, float *, int);
double covariance (float *, float *, int);
double rmse (float *, float *, int);
double bias (float *, float *, int);

/* matrix.c */
void freeDoubleVector(double * );
void freeFloatVector(float * );
void freeCharVector(unsigned char * );
void freeLongVector(long int * );
void freeShortVector(short int * );
double ** allocateDoubleMatrix(int, int);
float ** allocateFloatMatrix(int, int);
unsigned char ** allocateCharMatrix(int, int);
short int ** allocateShortMatrix(int, int);
long int ** allocateLongMatrix(int, int);
void freeDoubleMatrix(double **, int);
void freeFloatMatrix(float **, int);
void freeLongMatrix(long int **, int);
void freeShortMatrix(short int **, int);
void freeCharMatrix(unsigned char **, int);
double * allocateDoubleVector(int);
float * allocateFloatVector(int);
short int * allocateShortVector(int);
long int * allocateLongVector(int);
unsigned char * allocateCharVector(int);

double ** multiplyMatrix(double **, double **, int, int, int);
double ** trasposeMatrix(double **, int, int);
double * matrixToVector(double **, int, int);
double ** vectorToMatrix(double *, int, int);
double vectorLength(double *, int);
double ** addMatrix(double **, double **, int, int);
double ** subtractMatrix(double **, double **, int, int);
double ** mergeVector(double *, double *, int);
double ** createIdentityMatrix(int);
double ** scalarMultiplyMatrix(double **, double, int, int);
double ** diagonalMatrix(double **, int );
double * vectorInnerProduct(double *, double *, int);
double * vectorAddition(double *, double *, int);
double * vectorSubtract(double *, double *, int);
double * vectorDivision(double *, double *, int);
double * scalarMultiplyVector(double *, double, int);
double *  normalVector(double *, int);
double ** copyMatrix(double **, int, int);
double * copyVector(double *, int);
int printMatrix(double **, int, int);
int printVector(double *, int);
int resetMartixDoubleZero(double **, int, int);
int resetVectorDoubleZero(double *, int);

double matrixTrace(double **, int);
double * gaussianElimination (int, double **, double *);
double * backSubstitution (int, double **, double *);
double * forwardSubstitution (int, double **, double *);
double matrixDeterminant(double **, int);
double ** matrixInverse(double **, int);
double ** minorMatrix(double **, int, int, int);
double ** cofactorMatrix(double **, int);
int LUFactotisation (int, double **, double **, double **);
double * LUSolver(double **, double *, int);
double multipleCorrelation(double **, double *, int, int, double *);
double * multipleLinearRegression(double **, double *, int, int);
double * multiplyMatrixVector(double **, double *,  int, int);


/* root.c */
double secant(double (*f)(double), double, double);
double fixedPointIteration(double (*f)(double), double);
double bisect(double (*f)(double),double,double);
double newton(double (*f)(double), double (*fd)(double), double);
double regulaFalsi(double (*f)(double),double,double);
double brute (double (*f)(double), double, double, double);
double newtonSqaureRoot(double (*f)(double, double), double (*fd)(double), double, double);
double squareRoot(double);
double sqr(double, double);
double dsqr(double);

/* Fourier transform */
int complexConvert(float *, float *, float *, float *, int);
int RealDFT(float *, float *, int);
int ForwardDFT(float *, float *, float *, float *,  int);
int InverseDFT(float *, float *, float *, float *, int);

/* interpolation.c */
double linearInterpolationQuadratic(double *, double, double, double, double, double, double);
double linearInterpolation1D(double, double,  double, double, double);
double linearInterpolation(double *, double **, double *, int);
double ** allocateDoubleMatrix(int, int);
double * allocateDoubleVector(int);
double interp(int, int, double *, double *, double *,int);
double interp2D(int, double*, double *);
double mult(int, double, double *, double *);

/* integration.c */
double monteCarloIntegration1D(double (*f)(double), double, double, int);
double monteCarloIntegration(double (*f)(double *, int), double *, double *, int, int);
double trapizoidalRule(double (*f)(double), double, double, int);
double simpsonsRule(double (*f)(double), double, double, int);
double Integration(double (*f)(double), double (*method)(double (*f)(double), double, double, int), double, double);
double  trapizoidalRule2D(double (*f)(double, double), double, double, double, double, int, int);
double random01();
double random01Self(long int *r);

/* timeStep.c */
double euler(double (*f)(double, double), int, double);
double rungeKuttaFourthOrder(double (*f)(double, double), int, double);
double rungeKuttaSecondOrder(double (*f)(double, double), int, double);
double matsuno(double (*f)(double, double), int, double);
double leapFrog(double (*f)(double, double), int, double);

/* testStats.c */
double correlationTTest(double, int);
double rSignifcanceTest(double, int, int);
int correlationProbability(double, int);
double chiSquaredDistribution(double,int);
double chiProbability(double **, int, int);
double chiSquared(double **, int, int);

/* simplex.c */
int repositionSimplex(short int *, double **, int);
int sortSimplex(double *, int, short int *);
int calcCentroid(double **, int, double *);
double simplex (double (*f)(double *), double **, int);
double merit_function(double *, double *, double *, int, double (*f)(double *, double));
double simplex_fit (double **, int n, double *, double *, int, double (*f)(double *, double));


/* golden.c */
double golden (double (*f)(double), double, double);
double goldenFit (double, double, float *, float *, int, double (*f)(float, float));
float merit_function_golden(double *, float *, float *, int, double (*f)(float, float));

/* testFunctions.c */
double funcTimeStep(double, double);
double testFunctionSimplex(double *);
double exponent(double);
double cosine(double);
double poly(double);
double dpoly(double);
double functionTest1(double x);
double functionTest2(double x);
double functionTest3(double x);
double funcTest(double);
double funcTest5(double, double);
double funcTest6(double, double);
double testFunctionSimplexRosenbrock(double *);
double testFunctionSimplexFit(double *);
double test_linear_model(float, float);
double test_linear_model_2(double *p, double h);





