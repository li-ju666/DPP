#include "mat.h"

int main(int argc, char** argv){
    if(argc != 3){
	printf("Error: invalid parameters. Input file name and output file name required. \n"); 
	return -1; 
    }
    /* MPI initialization */
    MPI_status status; 
    int rank, size; 

    MPI_Init(int argc, char** argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

    /* Local variable declaration */
    int dim; 
    float *A_all, *B_all, *C_all, *A_local, *B_local, *A_cal, *B_cal, *C; 
    if(rank == 0){
	dim = read_input(argv[1], &A_all, &B_all); 
	if(dim <= 0){
	    printf("Error: Read input data failed! \n"); 
	    return -1; 
	}
	dim = dim/(int)sqrt((double)size); 
	MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    }
    
    A_local = malloc(dim*dim*sizeof(float)); 
    B_local = malloc(dim*dim*sizeof(float)); 

    if(A_local == NULL || B_local == NULL){
	printf("From rank %d: memory allocation for local storage matrix failed! \n", rank); 
	return -1; 
    }

    A_cal = malloc(dim*dim*sizeof(float)); 
    B_cal = malloc(dim*dim*sizeof(float)); 
    C = malloc(dim*dim*sizeof(float)); 

    if(A_cal == NULL || B_cal == NULL || C == NULL){
	printf("From rank %d: memory allocation for local calculation matrix failed! \n", rank); 
	return -1; 
    }

    

}

/* int main(int argc, char** argv){ */
/*     float *A, *B; */ 
/*     int dim; */ 
/*     dim = read_input(argv[1], &A, &B); */ 
/*     vis(A, dim); */ 
/*     vis(B, dim); */ 
/* } */
