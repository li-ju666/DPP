#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

#define ASSIGN_TAG 1111
#define SEND_TAG 1349
#define ARRAY_SEND_TAG 666

int cmpfunc(const void*, const void*); 
void vis(const int*, int); 
int cal_dim(int); 
int read_input(const char*, int**); 
int get_rank(const int*, int, int); 
int write_output(char*, const int*, int);

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
    MPI_Request request; 
    MPI_Status status;  
 
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */ 
    dim = cal_dim(size);

    int* range = malloc(dim*sizeof(int)); 
    int* period = malloc(dim*sizeof(int)); 
    int* my_coord = malloc(dim*sizeof(int)); 
    int* target_coord = malloc(dim*sizeof(int)); 
    int target_rank; 

    for(int i=0; i<dim; i++){
	range[i] = 2; 
	period[i] = 1; 
    } 
    
    MPI_Cart_create(MPI_COMM_WORLD, dim, range, period, 1, &CUBE_COMM);  
    MPI_Comm_rank(CUBE_COMM, &rank); 
    MPI_Cart_coords(CUBE_COMM, rank, dim, my_coord); 
    /* printf("From rank %d: my coord is %d, %d, %d, %d. \n", rank, */ 
	    /* my_coord[0], my_coord[1], my_coord[2], my_coord[3]); */ 

    /* Local variable declaration */
    int* array_all; 
    int* array_local; 
    int num_all, length;
    int* lengths, *orders;  

    if(rank == 0){
	num_all = read_input(argv[1], &array_all); 
	if(num_all < 1){
	    printf("Error: Invalid input file! \n"); 
	}
	/* vis(array_all, num_all); */ 
    }
    MPI_Bcast(&num_all, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    if(rank == size-1){
	length = num_all - num_all/size*(size-1); 
    }
    else{
	length = num_all/size; 
    }
    array_local = (int*)malloc(sizeof(int)*length); 
    if(array_local == NULL){
	printf("Error: Failed to allocate memory! \n"); 
    }

    if(rank == 0){
	memcpy(array_local, array_all, length*sizeof(int)); 
	for(int i=1; i<size; i++){
	    length = num_all/size; 
	    if(i == size-1){
		length = num_all - num_all/size*(size-1); 
	    }
	    /* printf("%d round: send length %d, index %d. ", i, length, num_all/size*i); */ 
	    MPI_Send(&(array_all[num_all/size*i]), length, MPI_INT, i, ASSIGN_TAG+i, MPI_COMM_WORLD); 
	}
	length = num_all/size; 
    }
    else{
	/* printf("Going to receive from rank %d ---- length %d. \n", rank, length); */ 
	MPI_Recv(array_local, length, MPI_INT, 0, ASSIGN_TAG+rank, MPI_COMM_WORLD, &status); 
	/* printf("Received from rank %d. \n", rank); */ 

    }

    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: ", rank); */ 
	    /* vis(array_local, length); */ 
	    /* printf("\n"); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 
    /* } */
    qsort(array_local, length, sizeof(int), cmpfunc); 



    int pivot, sum; 
    MPI_Comm SUB_COMM = CUBE_COMM; 
    int sub_rank = rank, sub_size = size; 
    int* medians, *buffer; 
    int smaller, largeq, buffer_size; 
    for(int i=0; i<dim; i++){
	memcpy(target_coord, my_coord, dim*sizeof(int)); 
	/* pivot choosing */
	switch(strategy){
	    case 1: 
		if(sub_rank == 0){
		    if(length > 0){
			pivot = array_local[length/2]; 
		    }
		    else{
			pivot = -1; 
		    }
		}
		MPI_Bcast(&pivot, 1, MPI_INT, 0, SUB_COMM); 
		break; 
	    case 2:
		/* pivot = ???; */
		if(sub_rank == 0){
		    medians = malloc(sizeof(int)*sub_size); 
		}
		MPI_Barrier(SUB_COMM); 
		MPI_Gather(&(array_local[length/2]), 1, MPI_INT, medians, 1, MPI_INT, 0, SUB_COMM); 
		if(sub_rank == 0){
		    /* vis(medians, sub_size); */ 
		    qsort(medians, sub_size, sizeof(int), cmpfunc);
		    /* vis(medians, sub_size); */  
		    pivot = medians[sub_size/2];
		    free(medians);  
		}
		MPI_Bcast(&pivot, 1, MPI_INT, 0, SUB_COMM); 
		break; 
	    case 3: 
		pivot = array_local[length/2]; 
		MPI_Allreduce(&pivot, &sum, 1, MPI_INT, MPI_SUM, SUB_COMM); 
		pivot = sum/sub_size; 
	}
    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: my pivot is %d. ", rank, pivot); */ 
	    /* printf("\n"); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 

    /* } */
	smaller = 0;
	for(int i=0; i<length; i++){
	    if(array_local[i] < pivot){
		smaller++; 
	    }
	    else{
		break; 
	    }
	}
	largeq = length - smaller; 
    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: my smaller is %d, larger equal is %d. ", rank, smaller, largeq); */ 
	    /* printf("\n"); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 
    /* } */
	memcpy(target_coord, my_coord, dim*sizeof(int)); 
	target_coord[i] = (target_coord[i]+1)%2; 
	MPI_Cart_rank(CUBE_COMM, target_coord, &target_rank); 
    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: my target rank is %d. \n ", rank, target_rank); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 
    /* } */
	if(my_coord[i]%2 == 0){
	    MPI_Isend(&largeq, 1, MPI_INT, target_rank, SEND_TAG+rank, CUBE_COMM, &request); 
	    MPI_Recv(&buffer_size, 1, MPI_INT, target_rank, SEND_TAG+target_rank, CUBE_COMM, &status); 
	    MPI_Wait(&request, &status);
	    buffer_size += smaller;  
	}
	else{
	    MPI_Isend(&smaller, 1, MPI_INT, target_rank, SEND_TAG+rank, CUBE_COMM, &request); 
	    MPI_Recv(&buffer_size, 1, MPI_INT, target_rank, SEND_TAG+target_rank, CUBE_COMM, &status); 
	    MPI_Wait(&request, &status); 
	    buffer_size += largeq; 
	}
    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: my buffer size is %d. ", rank, buffer_size); */ 
	    /* printf("\n"); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 

    /* } */

	buffer = malloc(sizeof(int)*buffer_size); 
	if(my_coord[i]%2 == 0){
	    memcpy(buffer, array_local, smaller*sizeof(int)); 
	    MPI_Isend(&(array_local[smaller]), largeq, MPI_INT, target_rank, ARRAY_SEND_TAG+rank, CUBE_COMM, &request); 
	    MPI_Recv(&(buffer[smaller]), buffer_size - smaller, MPI_INT, 
		    target_rank, ARRAY_SEND_TAG+target_rank, CUBE_COMM, &status); 
	    MPI_Wait(&request, &status); 
	}
	else{
	    memcpy(buffer, &(array_local[smaller]), largeq*sizeof(int)); 
	    MPI_Isend(array_local, smaller, MPI_INT, target_rank, ARRAY_SEND_TAG+rank, CUBE_COMM, &request); 
	    MPI_Recv(&(buffer[largeq]), buffer_size - largeq, MPI_INT, 
		    target_rank, ARRAY_SEND_TAG+target_rank, CUBE_COMM, &status); 
	    MPI_Wait(&request, &status); 
	}
	free(array_local); 
	array_local = buffer; 
	length = buffer_size; 
	qsort(array_local, length, sizeof(int), cmpfunc); 
    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: ", rank); */ 
	    /* vis(array_local, length); */ 
	    /* printf("\n"); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 

    /* } */
	/* do somgthing */
	MPI_Comm_split(SUB_COMM, my_coord[i], 0, &SUB_COMM); 
	MPI_Comm_rank(SUB_COMM, &sub_rank); 
	MPI_Comm_size(SUB_COMM, &sub_size); 
    }
    int order = 0; 
    for(int i=0; i<dim; i++){
	order += my_coord[i] << (dim-i-1); 
    }
    /* printf("From rank %d: my order is %d. \n", rank, order); */ 
    if(rank == 0){
	lengths = malloc(sizeof(int)*size); 
	orders = malloc(sizeof(int)*size); 
    }
    MPI_Gather(&length, 1, MPI_INT, lengths, 1, MPI_INT, 0, CUBE_COMM);
    MPI_Gather(&order, 1, MPI_INT, orders, 1, MPI_INT, 0, CUBE_COMM);
    MPI_Isend(array_local, length, MPI_INT, 0, ARRAY_SEND_TAG+rank, CUBE_COMM, &request); 
    if(rank == 0){
	int position = 0; 
	for(int i=0; i<size; i++){
	    target_rank = get_rank(orders, size, i); 
	    MPI_Recv(&(array_all[position]), lengths[target_rank], MPI_INT, 
		    target_rank, ARRAY_SEND_TAG+target_rank, CUBE_COMM, &status);
	    position += lengths[target_rank];  
	}
    }
    MPI_Wait(&request, &status); 
    if(rank == 0){
	/* vis(array_all, num_all); */
	printf("Running time to be measured! \n"); 
	write_output(argv[2], array_all, num_all); 	
    }
    /* for(int i=0; i<size; i++){ */
	/* if(rank == i){ */
	    /* printf("From rank %d: order - %d\n", rank, order); */ 
	    /* vis(array_local, length); */ 
	    /* printf("\n"); */ 
	/* } */
	/* MPI_Barrier(MPI_COMM_WORLD); */ 
    /* } */
    MPI_Finalize(); 
  
}

