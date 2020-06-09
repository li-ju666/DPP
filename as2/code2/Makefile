CC=mpicc

matmult: mat.c main.c mat.h
	$(CC) mat.c main.c -o matmult -O3 -std=c99 -lm 
clean: 
	rm -f mat.o main.o matmult
