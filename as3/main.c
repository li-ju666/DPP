#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

#define ASSIGN_TAG 1111
int cmpfunc(const void*, const void*); 
void vis(const int*, int); 
int cal_dim(int); 
int read_input(const char*, int**); 


int main(int argc, char** argv){
    /* Input parameters check */
    if(argc != 4){
	printf("Error: Invalid parameters. \nParameters: Input file Output file Strategy \n"); 
	return -1; 
    }
    
    int strategy = atoi(argv[3]); 
    if(strategy < 0 || strategy > 3){
	printf("Error: Invalid pivot selection strategy! Please choose from 1, 2 and 3. \n"); 
    }

    /* MPI initialization */
    int rank, size, dim; 
    MPI_Comm CUBE_COMM; 
 
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */ 
    dim = cal_dim(size);

    int* range = malloc(dim*sizeof(int)); 
    int* period = malloc(dim*sizeof(int)); 
    int* my_coord = malloc(dim*sizeof(int)); 
    int* target_coord = malloc(dim*sizeof(int)); 

    for(int i=0; i<dim; i++){
	range[i] = 2; 
	period[i] = 1; 
    } 
    
    MPI_Cart_create(MPI_COMM_WORLD, dim, range, period, 1, &CUBE_COMM);  
    MPI_Comm_rank(CUBE_COMM, &rank); 
    MPI_Cart_coords(CUBE_COMM, rank, dim, my_coord); 
    printf("From rank %d: my coord is %d, %d, %d, %d. \n", rank, 
	    my_coord[0], my_coord[1], my_coord[2], my_coord[3]); 

    /* Local variable declaration */
    int total_num; 
    int* array_all; 
    int* array_local; 
    int num_all; 
    MPI_Status status; 

    if(rank == 0){
	num_all = read_input(argv[1], &array_all); 
	if(num_all < 1){
	    printf("Error: Invalid input file! \n"); 
	}
	vis(array_all, num_all); 
    }
    MPI_Bcast(&num_all, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    array_local = (int*)malloc(2*sizeof(int)*num_all/size); 
    if(array_local == NULL){
	printf("Error: Failed to allocate memory! \n"); 
    }

    if(rank == 0){
	memcpy(array_local, array_all, num_all/size*sizeof(int)); 
	for(int i=1; i<size; i++){
	MPI_Send(&(array_all[i*num_all/size]), num_all/size, MPI_INT, i, ASSIGN_TAG+i, MPI_COMM_WORLD); 
	}
    }
    else{
	MPI_Recv(array_local, num_all/size, MPI_INT, 0, ASSIGN_TAG+rank, MPI_COMM_WORLD, &status); 
    }

qsort(array_local, num_all/size, sizeof(int), cmpfunc); 

    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: ", rank); */ 
	    /* vis(array_local, num_all/size); */ 
	    /* printf("\n"); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 

    /* } */
    int pivot; 



    MPI_Finalize(); 
    
}

