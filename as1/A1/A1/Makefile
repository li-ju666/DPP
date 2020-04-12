CC = mpicc
FLAGS = -O3
STD = -std=c99

all: mpi_main.c stencil.h
	$(CC) mpi_main.c $(FLAGS) $(STD) -o stencil

clean:  
	rm -f stencil
