#####################################
#  Makefile for compiling Numerical.
#  Creates both a static and a
#  dynamic library. 
#
# Need to set up LD_LIBRARY_PATH
# for shared libraries.
#
# Will Grey
# 14 June 2015
#####################################

CC=gcc
BINDIR = bin/
SRCDIR = src/
LIBDIR = lib/
LIB = -lm
CFLAGS = -Wall -I.

SRC = matrix.c univariateStats.c utility.c sequence.c\
      fourierTransform.c interpolation.c integration.c root.c\
      timeStep.c statsTest.c simplex.c testFunctions.c bivariateStats.c golden.c

CSources = $(addprefix $(SRCDIR),$(SRC))

all : $(BINDIR)numerical

# Use numerical.c to test the library. 
$(BINDIR)numerical : $(LIBDIR)libnumericalstatic.a $(LIBDIR)libnumerical.so numerical.c
	$(CC) numerical.c -o $(BINDIR)numerical -lnumericalstatic -Llib $(LIB) $(CFLAGS)

# Create a static library
$(LIBDIR)libnumericalstatic.a : $(CSources) 
	$(CC) $(CSources) -c $(LIB) $(CFLAGS)
	ar -cq $(LIBDIR)libnumericalstatic.a *.o
	rm *.o

# Create a shared library
$(LIBDIR)libnumerical.so : $(CSources) 
	$(CC) $(CSources) -c $(LIB) -fPIC $(CFLAGS)
	$(CC) -shared -o $(LIBDIR)libnumerical.so *.o
	rm *.o

#$(BINDIR)numerical : $(CSources) 
#	$(CC) $(CSources)\
#       -o $(BINDIR)numerical $(CFLAGS) $(LIB)

clean :
	rm -f $(BINDIR)numerical
	rm -f $(LIBDIR)libnumericalstatic.a
	rm -f $(LIBDIR)libnumerical.so
	rm -f *.o
        
