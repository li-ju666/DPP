###############################################################################
# Makefile for assignment 1, Parallel and Distributed Computing 2020.
###############################################################################

CC = mpicc
# CFLAGS = -std=c99 -g -O3
CFLAGS = -g -O3
LIBS = -lm

BIN = stencil

all: $(BIN)

stencil: stencil.c stencil.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN)
