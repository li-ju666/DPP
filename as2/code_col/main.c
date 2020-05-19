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

    int grid_dim[2], grid_period[2], my_coord[2] = {0,0}, 
	up_coord[2], down_coord[2], 
	left, right, up, down; 
    char* output_file = argv[2];
    /* MPI initialization */
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    MPI_Request Arequest[size], Brequest[size]; 
    MPI_Status status[size]; 

    grid_dim[0] = grid_dim[1] = (int)sqrt(size); 
    grid_period[0] = grid_period[1] = 1; 
    
    /* if(rank == 0){printf("The grid dim are %d, %d. \n", grid_dim[0], grid_dim[1]); } */
    if(grid_dim[0]*grid_dim[1] != size){
	printf("Error: the number of cores cannot be grided. \n"); 
	MPI_Abort(MPI_COMM_WORLD, 1); 
	return -1; 
    }
    
    /* Create topology */
    MPI_Cart_create(MPI_COMM_WORLD, 2, grid_dim, grid_period, 1, &GRID_COMM); 
    MPI_Cart_coords(GRID_COMM, rank, grid_dim[0], my_coord);
    
    /* Find neighbours */

    up_coord[1] = down_coord[1] = my_coord[1]; 
    up_coord[0] = (my_coord[0] - 1 + grid_dim[0]) % grid_dim[0]; 
    down_coord[0] = (my_coord[0] + 1) % grid_dim[0]; 

    MPI_Cart_rank(GRID_COMM, up_coord, &up); 
    MPI_Cart_rank(GRID_COMM, down_coord, &down);

    /* Local variable declaration */
    int dim_all, dim; 
    float *A_all, *B_all, *C_all, *A_local, *A_cal, *B_cal, *C; 
    if(rank == 0){
	dim_all = read_input(argv[1], &A_all, &B_all); 
	if(dim_all <= 0){
	    printf("Error: Read input data failed! \n"); 
	    return -1; 
	}
	C_all = malloc(sizeof(float)*dim_all*dim_all); 
	/* printf("Original data dim is: %d. \n", dim); */
	dim = dim_all/(int)sqrt(size); 
	if(dim * (int)sqrt(size) != dim_all){
	    printf("The matrix cannot be grided. \n"); 
	    MPI_Abort(GRID_COMM, 1); 
	    return -1; 
	}
    }
    
    MPI_Bcast(&dim, 1, MPI_INT, 0, GRID_COMM); 
    /* printf("From process %d: local data dim is: %d. \n", rank, dim); */   
    
    /* Define checkboard data type (for task assigning) */
    /* MPI_Type_vector(dim, dim,dim*grid_dim[0], MPI_FLOAT, &checkboard_type); */ 
    /* MPI_Type_commit(&checkboard_type); */ 
    
    /* Another approach - cyclic checkboard split */
    int block_len[dim]; 
    MPI_Aint indices[dim]; 
    MPI_Datatype data_type[dim];
     
    for(int i=0; i<dim; i++){
	block_len[i] = 1; 
	indices[i] = i*grid_dim[0]*sizeof(float);
	data_type[i] = MPI_FLOAT; 
    }
    /* indices[dim-1] = grid_dim[0]*dim*grid_dim[0]*sizeof(float); */ 
    /* if(rank == 0){ */
	/* printf("Dim is %d. \n===========================\n", dim); */ 
	/* for(int i=0; i<dim; i++){ */
	    /* printf("Index is %d. \n", indices[i]); */ 
	/* } */
    /* } */
    MPI_Type_create_struct(dim, block_len, indices, data_type, &checkboard_type); 
    MPI_Type_create_resized(checkboard_type, 0, grid_dim[0]*grid_dim[0]*dim*sizeof(float), &checkboard_type);
    MPI_Type_commit(&checkboard_type); 
    

    /* Allocate memory for local stored partial-A and calculation-required partial A and B */
    A_local = malloc(dim*dim*sizeof(float)); 
    A_cal = malloc(dim*dim*sizeof(float)); 
    B_cal = malloc(dim*dim*sizeof(float)); 
    C = calloc(dim*dim, sizeof(float)); 

    if(A_local == NULL || A_cal == NULL || B_cal == NULL || C == NULL){
	printf("From rank %d: memory allocation for matrices failed! \n", rank); 
	return -1; 
    }
    
    /* Assign matrix A and B to all processes */ 
    if(rank == 0){
	for(int i=0; i<size; i++){
	    /* printf("Index is: %d. \n", i%grid_dim[0]*dim*grid_dim[0] + i/grid_dim[0]); */ 
	    MPI_Isend(&(A_all[(i%grid_dim[0])*dim*grid_dim[0] + i/grid_dim[0]]),  
		    dim, checkboard_type, i, A_TAG+i, GRID_COMM, &Arequest[i]); 
	    MPI_Isend(&(B_all[(i%grid_dim[0])*dim*grid_dim[0] + i/grid_dim[1]]),  
		    dim, checkboard_type, i, B_TAG+i, GRID_COMM, &Brequest[i]); 
	}
    }
    MPI_Recv(A_local, dim*dim, MPI_FLOAT, 0, A_TAG+rank, GRID_COMM, &status[0]); 
    MPI_Recv(B_cal, dim*dim, MPI_FLOAT, 0, B_TAG+rank, GRID_COMM, &status[0]); 
    if(rank == 0){
	MPI_Waitall(size, Arequest, status); 
	MPI_Waitall(size, Brequest, status); 
    }
    
    /* Define Row communicator */
    MPI_Comm ROW_COMM;  
    MPI_Comm_split(GRID_COMM, my_coord[0], rank, &ROW_COMM); 
    
    /* Calculation for A*B parallel*/
    for(int i=0; i<grid_dim[0]; i++){
	/* Broadcast target partial matrix A */
	if((my_coord[0]+i)%grid_dim[0] == my_coord[1]){
	    memcpy(A_cal, A_local, sizeof(float)*dim*dim); 
	}
	MPI_Bcast(A_cal, dim*dim, MPI_FLOAT, 
		(my_coord[0]+i)%grid_dim[0], ROW_COMM);
	/* Local matrix calculation */
	float* temp = multiply(A_cal, B_cal, dim); 
	sum(C, temp, dim);
	free(temp);
	/* vis(A_local, dim); */ 
	/* printf("===========================\n"); */ 	
	/* vis(A_cal, dim); */ 
	/* printf("====================\n"); */
	/* vis(B_cal, dim); */ 
	/* Partial matrix B scrolling up */
	MPI_Isend(B_cal, dim*dim, MPI_FLOAT, up, B_TAG+rank, GRID_COMM, &Brequest[0]); 
	MPI_Recv(B_cal, dim*dim, MPI_FLOAT, down, B_TAG+down, GRID_COMM, &status[0]); 
	MPI_Wait(&Brequest[0], &status[0]);
	MPI_Barrier(GRID_COMM); 	
    }
    /* Collect result from all processes to process with rank 0 */ 
    MPI_Isend(C, dim*dim, MPI_FLOAT, 0, B_TAG+rank, GRID_COMM, &Brequest[0]); 
    if(rank == 0){
	for(int i=0; i<size; i++){
	    MPI_Recv(&(C_all[(i%grid_dim[0])*dim*grid_dim[0] + i/grid_dim[0]]), 
		    dim, checkboard_type, i, B_TAG+i, GRID_COMM, &status[0]); 
	}
	/* Print result out */
	printf("Matrix A: \n"); 
	vis(A_all, dim*grid_dim[0]); 
	printf("Matrix B: \n"); 
	vis(B_all, dim*grid_dim[0]); 
	printf("A * B: \n"); 
	vis(C_all, dim*grid_dim[0]);
	if(write_output(output_file, C_all, dim*grid_dim[0]) != 0){
	    printf("Write output matrix failed! \n"); 
	}
	free(C_all); 
	free(A_all); 
	free(B_all); 
    }
    MPI_Wait(&Brequest[0], &status[0]); 
    free(A_local); 
    free(A_cal); 
    free(B_cal); 
    free(C); 
    MPI_Type_free(&checkboard_type); 
    MPI_Finalize();
    return 0;  
}
