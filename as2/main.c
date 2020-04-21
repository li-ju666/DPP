#include "mat.h"
#include <math.h>
#include <mpi.h>
#include <string.h>
#define A_TAG 1349
#define B_TAG 204

int main(int argc, char** argv){
    if(argc != 3){
	printf("Error: invalid parameters. Input file name and output file name required. \n"); 
	return -1; 
    }
    /* MPI initialization */
    int rank, size; 
    MPI_Comm GRID_COMM; 
    MPI_Datatype checkboard_type;

    int grid_dim[2], grid_period[2], my_coord[2], 
	left_coord[2], right_coord[2], 
	up_coord[2], down_coord[2], 
	left, right, up, down; 
    
    /* MPI initialization */
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    MPI_Request Arequest[size], Brequest[size]; 
    MPI_Status status[size]; 

    grid_dim[0] = grid_dim[1] = (int)sqrt(size); 
    grid_period[0] = grid_period[1] = 1; 

    if(grid_dim[0]*grid_dim[1] != size){
	printf("Error: the number of cores cannot be grided. \n"); 
	MPI_Abort(MPI_COMM_WORLD, 1); 
	return -1; 
    }
    
    /* Create topology */
    MPI_Cart_create(MPI_COMM_WORLD, 2, grid_dim, grid_period, 1, &GRID_COMM); 
    MPI_Cart_coords(GRID_COMM, rank, grid_dim[0], my_coord);
    
    /* Find neighbours */
    left_coord[0] = right_coord[0] = my_coord[0]; 
    left_coord[1] = (my_coord[1] - 1 + grid_dim[1]) % grid_dim[1]; 
    right_coord[1] = (my_coord[1] + 1) % grid_dim[1]; 

    up_coord[1] = down_coord[1] = my_coord[1]; 
    up_coord[0] = (my_coord[0] - 1 + grid_dim[0]) % grid_dim[0]; 
    down_coord[0] = (my_coord[0] + 1) % grid_dim[0]; 

    MPI_Cart_rank(GRID_COMM, left_coord, &left); 
    MPI_Cart_rank(GRID_COMM, right_coord, &right); 
    MPI_Cart_rank(GRID_COMM, up_coord, &up); 
    MPI_Cart_rank(GRID_COMM, down_coord, &down);

    /* printf("From process %d: left: %d, right: %d, up: %d, down: %d. \n", */ 
	    /* rank, left, right, up, down); */ 

    /* Local variable declaration */
    int dim; 
    float *A_all, *B_all, *C_all, *A_local, *B_local, *A_cal, *B_cal, *C; 
    if(rank == 0){
	dim = read_input(argv[1], &A_all, &B_all); 
	if(dim <= 0){
	    printf("Error: Read input data failed! \n"); 
	    return -1; 
	}
	/* printf("Original data dim is: %d. \n", dim); */ 
	dim = dim/(int)sqrt(size); 
    }
    
    MPI_Bcast(&dim, 1, MPI_INT, 0, GRID_COMM); 
    /* printf("From process %d: local data dim is: %d. \n", rank, dim); */   
    
    /* Define checkboard data type (for task assigning) */
    MPI_Type_vector(dim*dim, 1, size, MPI_FLOAT, &checkboard_type); 
    MPI_Type_commit(&checkboard_type); 
    
    /* Allocate memory for local storage data and calculation-required data */
    A_local = malloc(dim*dim*sizeof(float)); 
    B_local = malloc(dim*dim*sizeof(float)); 

    if(A_local == NULL || B_local == NULL){
	printf("From rank %d: memory allocation for local storage matrix failed! \n", rank); 
	return -1; 
    }

    A_cal = malloc(dim*dim*sizeof(float)); 
    B_cal = malloc(dim*dim*sizeof(float)); 
    C = calloc(dim*dim, sizeof(float)); 

    if(A_cal == NULL || B_cal == NULL || C == NULL){
	printf("From rank %d: memory allocation for local calculation matrix failed! \n", rank); 
	return -1; 
    }
    
    if(rank == 0){
	for(int i=0; i<size; i++){
	    MPI_Isend(&(A_all[i]), 1, checkboard_type, i, A_TAG+i, GRID_COMM, &Arequest[i]); 
	    MPI_Isend(&(B_all[i]), 1, checkboard_type, i, B_TAG+i, GRID_COMM, &Brequest[i]); 
	}
    }
    MPI_Recv(A_local, dim*dim, MPI_FLOAT, 0, A_TAG+rank, GRID_COMM, &status[0]); 
    MPI_Recv(B_local, dim*dim, MPI_FLOAT, 0, B_TAG+rank, GRID_COMM, &status[0]); 
    if(rank == 0){
	MPI_Waitall(size, Arequest, status); 
	MPI_Waitall(size, Brequest, status); 
    }
    /* printf("From process %d: \n", rank); */  
    /* vis(A_local, dim); */ 
    /* vis(B_local, dim); */ 
    MPI_Comm ROW_COMM, COL_COMM; 
    MPI_Comm_split(GRID_COMM, my_coord[0], rank, &ROW_COMM); 
    MPI_Comm_split(GRID_COMM, my_coord[1], rank, &COL_COMM); 
    
    memcpy(B_cal, B_local, sizeof(float)*dim*dim); 

    for(int i=0; i<dim; i++){
	if(my_coord[0]+i == my_coord[1]){
	    memcpy(A_cal, A_local, sizeof(float)*dim*dim); 
	}
	MPI_Bcast(A_cal, sizeof(float)*dim*dim, MPI_FLOAT, 
		my_coord[0]*dim+(my_coord[0]+i)%dim, ROW_COMM);
	sum(C, multiply(A_cal, B_cal, dim), dim); 
	MPI_Isend(B_local, dim*dim, MPI_FLOAT, up, B_TAG+rank, GRID_COMM, &Brequest[0]); 
	MPI_Recv(B_cal, dim*dim, MPI_FLOAT, down, B_TAG+down, GRID_COMM, &status[0]); 
	MPI_Wait(&Brequest[0], &status[0]); 
    }
    
    printf("From process %d: ", rank); 
    vis(C, dim); 

    MPI_Finalize();
    return 0;  
}


