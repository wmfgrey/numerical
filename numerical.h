
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>

#define MAXIMUM(x,y) ((x) > (y) ? (x) : (y) )
#define MINIMUM(x,y) ((x) < (y) ? (x) : (y) )
#define SQUARE(x) ((x) * (x))
#define ARRAY_LENGTH(x)  (sizeof(x) / sizeof((x)[0]))
#define ABS(x) ((x) >= 0 ? (x) : (-x))

#define PI 3.1415926535897932
#define E 2.7182818284590452
#define MAX_PATH_NAME_LENGTH 1000 
#define NUM_DISTRIBUTION_BINS 100

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
double monteCarloIntegrationPi();
unsigned char *prime(int);
long *fibonacci(int);
double bayes(double, double, double);
pdf* normalDistribution(double, double, int, double, int);
double gammaFunction(double k);
double chiSquaredDistribution(double, int);

/* bivariateStats.c */
regressionCoefficients regression(float *, float *, int);
double correlation (float *, float *, int);
double covariance (float *, float *, int);
double rmse (float *, float *, int);
double bias (float *, float *, int);

/* testStats.c */

/* matrix.c */

/* root.c */




	
