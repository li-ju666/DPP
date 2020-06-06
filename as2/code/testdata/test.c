#include <stdio.h>
#include "mat.h"
#include <mpi.h>

int main(int argc, char** argv){
	int dim=3; 
	/* float *A, *B, *C; */ 
	/* dim = read_input(argv[1], &A, &B); */ 
	/* double start; */ 
	/* start = MPI_Wtime(); */ 
	/* C = multiply(A, B, dim); */
	/* printf("%f\n", MPI_Wtime()-start); */ 
	/* free(A); */ 
	/* free(B); */ 
	/* free(C); */
	float A[9] = {1,2,3,4,5,6,7,8,9}; 
	float B[9] = {2,3,4,5,6,7,8,9,10}; 
	vis(A, 3); 
	vis(B, 3); 
	float* C = multiply(A, B, dim);
	vis(C, 3); 
	free(C); 

}
