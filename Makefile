#####################################
# Produces both static and shared 
# libraries for numerical routines.
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
      timeStep.c statsTest.c simplex.c testFunctions.c bivariateStats.c

CSources = $(addprefix $(SRCDIR),$(SRC))

all : $(BINDIR)numerical

$(BINDIR)numerical : $(LIBDIR)libnumericalstatic.a $(LIBDIR)libnumerical.so numerical.c
	$(CC) numerical.c -o $(BINDIR)numerical -lnumerical -Llib $(LIB) $(CFLAGS)

$(LIBDIR)libnumericalstatic.a : $(CSources) 
	$(CC) $(CSources) -c $(LIB) $(CFLAGS)
	ar -cq $(LIBDIR)libnumericalstatic.a *.o
	rm *.o

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
        




