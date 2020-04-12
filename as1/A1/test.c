#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[]){
    
    MPI_Init(NULL, NULL);
    int size, rank; 
    
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    printf("Me myself is: %d of %d. \n", rank, size);
    
    MPI_Finalize();  
}
