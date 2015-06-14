#####################################
# Makefile for compiling Numerical.
# Will Grey
# 14 June 2015
#####################################

CC=gcc
BINDIR = bin/
SRCDIR = src/
LIB = -lm
CFLAGS = -Wall -I.

SRC=numerical.c univariateStats.c utility.c sequence.c bivariateStats.c
CSources=$(addprefix $(SRCDIR),$(SRC))

all : $(BINDIR)numerical

$(BINDIR)numerical : $(CSources) 
	$(CC) $(CSources)\
       -o $(BINDIR)numerical $(CFLAGS) $(LIB)

clean :
	rm -f $(BINDIR)numerical




