CC=mpicc

all: matmul

matmul: mat.c main.c mat.h
	$(CC) mat.c main.c -o matmul -O3 -std=c99 -lm

clean: 
	rm -f mat.o main.o matmul
