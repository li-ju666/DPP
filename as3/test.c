#include <stdio.h>
#include <mpi.h>

int main(){
    int rank, size, sub_rank; 
    MPI_Comm SUB_COMM; 

    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    SUB_COMM = MPI_COMM_WORLD; 
    sub_rank = rank; 
    for(int i=0; i<2; i++){
	MPI_Comm_split(SUB_COMM, sub_rank%2, 0, &SUB_COMM); 
	MPI_Comm_rank(SUB_COMM, &sub_rank); 
	printf("Loop %d --- From rank %d: my subrank is %d. \n", i, rank, sub_rank); 
    }
    MPI_Finalize(); 
}
