
/*
NAME:		utility.c 
DESCRIPTION: 	Utility tools for memory management,
                reading and writing to files and other
                housekeping.
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

FILE *openFile(char *path, char *fileName, char *rw)
{
 
 FILE *fp;
 char fullPathName[MAX_PATH_NAME_LENGTH];
 
 strcpy(fullPathName,path);
 strcat(fullPathName,fileName); 
 
 
 if ( !strcmp(rw, "rb") ){
 
  fp = fopen(fullPathName, "rb");
  if ( fp == NULL ){
   fprintf( stderr, "Error opening file %s\n", fullPathName);
   exit(EXIT_FAILURE); 
  }
  
 }
 
 if ( !strcmp(rw, "wb") ){
 
  fp = fopen(fullPathName, "wb");
  if ( fp == NULL ){
   fprintf( stderr, "Error creating %s\n", fullPathName);
   exit(EXIT_FAILURE); 
  }
 
 }
 
 if ( !strcmp(rw, "r") ){
 
  fp = fopen(fullPathName, "r");
  if ( fp == NULL ){
   fprintf( stderr, "Error opening file %s\n", fullPathName);
   exit(EXIT_FAILURE); 
  }
  
 }
 
 if ( !strcmp(rw, "w") ){
 
  fp = fopen(fullPathName, "w");
  if ( fp == NULL ){
   fprintf( stderr, "Error creating %s\n", fullPathName);
   exit(EXIT_FAILURE); 
  }
 
 }

 if ( !strcmp(rw, "a") ){
 
  fp = fopen(fullPathName, "a");
  if ( fp == NULL ){
   fprintf( stderr, "Error creating %s\n", fullPathName);
   exit(EXIT_FAILURE); 
  }
 
 }

 return fp;

}

float *readFloatData(char *pathName, char *fileName, int n)
{
 
 FILE *fp;
 float *data;

 fp=openFile(pathName, fileName, "rb");
 
 if((data = (float *) calloc(n,sizeof(float)))==NULL) memoryCheck();  
 fread(data, sizeof(float), n, fp);
 
 fclose(fp);
 
 return data;

} 

int memoryCheck(void)	
{	
 fprintf(stderr,"\nUnable to allocate memory.\n\n");
 exit(EXIT_FAILURE);
}

int createEnviHeader(char *pathName, char *fileName, int xdim, int ydim, int bands, int dataType){
 
 FILE *fp;
 fp=openFile(pathName, fileName, "w");
 
 fprintf(fp,"ENVI\n");
 fprintf(fp,"samples   = %d\n", xdim);
 fprintf(fp,"lines     = %d\n", ydim);
 fprintf(fp,"bands     = %d\n", bands);
 fprintf(fp,"data type = %d\n", dataType);
 fprintf(fp,"byte order = 1\n");
 
 fclose(fp);
 
 return 0;

}

int counter(unsigned int element, unsigned int nElements)

/* In code use something like:
 if (element % nElements == 0) counter (element,nElements);
*/

{
 unsigned int percentage;
	
 percentage = (element*100)/nElements;
 fprintf(stderr,"\rprocessing line: %d (%3d%% completed)",element, percentage);
 fflush(NULL);
		
 if(percentage == 100){
  fprintf(stderr,"\n"); 
  fflush(NULL);
 }	
	
 return(EXIT_SUCCESS);	 

}

double timer()
{
/* in code write:
 double startTime=timer();
at star of program

At end of program write:
 double endTime=timer();
 fprintf(stderr," %s: Runtime %f seconds\n", argv[0], endTime-startTime);
*/

 struct timeval tv;
 struct timezone tz;
 double walltime;

  gettimeofday(&tv, &tz);
  walltime = tv.tv_sec +tv.tv_usec*1.0e-6;

  return(walltime);
}


float * convertByte2Float(unsigned char *dataIn, int n){

 int i;
 float *dataOut;

 if((dataOut  = (float *) calloc(n, sizeof(float)))== NULL) memoryCheck();

 for (i=0;i<n;i++)
   *(dataOut +i) = (float)*(dataIn +i); 

 return dataOut;

}

float * convertDouble2Float(double *dataIn, int n){

 int i;
 float *dataOut;

 if((dataOut  = (float *) calloc(n, sizeof(float)))== NULL) memoryCheck();

 for (i=0;i<n;i++)
   *(dataOut +i) = (float)*(dataIn +i); 

 return dataOut;

}

float * convertLong2Float(unsigned long int *dataIn, int n){

 int i;
 float *dataOut;

 if((dataOut  = (float *) calloc(n, sizeof(float)))== NULL) memoryCheck();

 for (i=0;i<n;i++)
   *(dataOut +i) = (float)*(dataIn +i); 

 return dataOut;

}

float * convertShort2Float(unsigned short int *dataIn, int n){

 int i;
 float *dataOut;
 

 if((dataOut  = (float *) calloc(n, sizeof(float)))== NULL) memoryCheck();

 for (i=0;i<n;i++)
   *(dataOut +i) = (float)*(dataIn +i); 

 return dataOut;

}


int byteSwap(void *dataIn, int n,  int numBytes){
 
 /* convert from little endian to big endian or vice versa */
 /* 2 bytes: bytes ab reordered to ba */
 /* 4 bytes: bytes abcd reordered to dcba */ 
 /* 8 bytes: bytes abcdefgh reordered to hgfedcba */ 

 int i,j;
 unsigned char *value;
 unsigned char byte;
 int swap;
 swap = numBytes / 2;

 for (i=0;i<n*numBytes;i=i+numBytes){

  value = (unsigned char *)(dataIn+i); 
  
  for (j=0; j<swap;j++){
    
   byte = *(value + j);
   *(value) = *(value+(numBytes-1)-j);
   *(value+(numBytes-1)-j) = byte;
   
  }
 }

 return (EXIT_SUCCESS);

}



