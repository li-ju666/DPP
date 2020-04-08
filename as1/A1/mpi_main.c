#include "stencil.h"
#include <string.h>

#define LEFTTAG 1
#define RIGHTTAG 1111

int main(int argc, char* argv[]){
    if(argc != 4){
	printf("Usage: stencil input_file output_file number_of_applications\n");
	return 1;
    }
    char* input_name = argv[1]; 
    char* output_name = argv[2]; 
    int num_steps = atoi(argv[3]); 

    /* MPI initialization */
    MPI_Status status; 
    int rank, size, left, right, num_values; 
    double* localin, *localout, *alldata; 
    
    MPI_Init(NULL, NULL); 
    
    /* Assign local variables */
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    left = (rank-1+size)%size; 
    right = (rank+1+size)%size; 

    double leftbuff[2], rightbuff[2]; 
    MPI_Request lsendr, rsendr, lrecvr, rrecvr; 

    /* Read input file from Process 0*/
    if(rank == 0){
	/* double* input; */ 
	num_values = read_input(input_name, &alldata); 
	if(num_values < 0){return 2; }
    }

    /* Broadcast the number of values: to make other processes allocate memory */
    MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    /* printf("From process %d: My left is %d, my right is %d. \n", rank, left, right); */ 

    double h = 2.0*PI/num_values;
    const int STENCIL_WIDTH = 5;
    const int EXTENT = STENCIL_WIDTH/2;
    const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};
    printf("From process %d -- Stencil are: %.3f, %.3f, %.3f, %.3f, %.3f. \n\n", 
		rank, STENCIL[0], STENCIL[1], STENCIL[2], STENCIL[3], STENCIL[4]); 

    /* Allocate memory for local data */
    localin = (double*)malloc(num_values/size * sizeof(double)); 
    localout = (double*)malloc(num_values/size * sizeof(double)); 
    if(localin == NULL || localout == NULL){
	printf("Process %d: Could not allocate memory! \n", rank); 
    }

    /* Assign data (by p0) */
    if(rank == 0){
	/* Copy p0's own data */
	memcpy(localin, alldata, num_values/size*sizeof(double)); 
	/* Send data to other proceseses */
	for(int i=1; i<size; i++){
	    MPI_Ssend(&alldata[i*num_values/size], 
		num_values/size, MPI_DOUBLE, i, 111*i, MPI_COMM_WORLD); 
	}
    }
    /* Other process wait for assigned data */
    else{
	MPI_Recv(localin, num_values/size, MPI_DOUBLE, 0, 111*rank, MPI_COMM_WORLD, &status); 
    }
    
    for(int step=0; step<num_steps; step++){
	/* Communicate with neighbors for data */
	MPI_Isend(localin, EXTENT, MPI_DOUBLE, left, LEFTTAG*rank, MPI_COMM_WORLD, &lsendr); 
	MPI_Isend(&localin[num_values/size-2], EXTENT, MPI_DOUBLE, right, RIGHTTAG*rank, MPI_COMM_WORLD, &rsendr); 
	MPI_Irecv(leftbuff, EXTENT, MPI_DOUBLE, left, RIGHTTAG*left, MPI_COMM_WORLD, &lrecvr); 
	MPI_Irecv(rightbuff, EXTENT, MPI_DOUBLE, right, LEFTTAG*right, MPI_COMM_WORLD, &rrecvr); 

	MPI_Wait(&lsendr, &status); 
	MPI_Wait(&rsendr, &status); 
	MPI_Wait(&lrecvr, &status); 
	MPI_Wait(&rrecvr, &status);

	/* printf("From process %d: localin: %.3f, %.3f, %.3f, %.3f. \nLeft buff: %.3f, %.3f. \nRight buff: %.3f, %.3f. \n\n", */ 
	/* 	rank, localin[0], localin[1], localin[2], localin[3], */ 
	/* 	leftbuff[0], leftbuff[1], rightbuff[0], rightbuff[1]); */ 

	/* Apply stencil on local data */
	for(int i=0; i<EXTENT; i++){
	    double result=0; 
	    for(int j=0; j<STENCIL_WIDTH; j++){
		double base = 0; 
		if(i+j-EXTENT<0){ base = leftbuff[i+j]; }
		else{ base = localin[i+j-EXTENT]; }
		result += STENCIL[j] * base;
	    }
	    localout[i] = result; 
	}

	for(int i=EXTENT; i<num_values/size-EXTENT; i++){
	    double result = 0; 
	    for(int j=0; j<STENCIL_WIDTH; j++){
		result += STENCIL[j] * localin[i+j-EXTENT]; 
	    }
	    localout[i] = result; 
	}

	for(int i=num_values/size-EXTENT; i<num_values/size; i++){
	    double result = 0; 
	    for(int j=0; j<STENCIL_WIDTH; j++){
		double base = 0; 
		if(i+j-EXTENT>=num_values/size){ base = rightbuff[i+j-EXTENT-num_values/size]; }
		else{ base = localin[i+j-EXTENT]; }
		result += STENCIL[j] * base; 
	    }
	    localout[i] = result; 
	}
	/* printf("From process %d: localout: %.3f, %.3f, %.3f, %.3f. \n", */ 
	/* 	rank, localout[0], localout[1], localout[2], localout[3]); */ 
	MPI_Barrier(MPI_COMM_WORLD);

	/* Swap input and output */
	if(step < num_steps-1){
	   double* temp = localin; 
	   localin = localout; 
	   localout = temp; 
	}
	else{
	    printf("From process %d: result: %.3f, %.3f, %.3f, %.3f. \n", 
		rank, localout[0], localout[1], localout[2], localout[3]); 
	    
	    if(rank == 0){
		memcpy(alldata, localout, num_values/size*sizeof(double)); 
		for(int i=1; i<size; i++){
		    MPI_Recv(&alldata[i*num_values/size], num_values/size, 
			    MPI_DOUBLE, i, 111*i, MPI_COMM_WORLD, &status);
		   printf("Received data from %d: %.3f, %.3f, %.3f, %.3f. \n", 
			   alldata[i*num_values/size], alldata[i*num_values/size+1],alldata[i*num_values/size+2],alldata[i*num_values/size+3],i);  
		}
	    }
	    else{
		MPI_Ssend(localout, num_values/size, MPI_DOUBLE, 0, 111*rank, MPI_COMM_WORLD); 
	    }
	}
	/* MPI_Barrier(MPI_COMM_WORLD); */
    }
    /* printf("From process %d: left buff: %.2f, %.2f. right buff: %.2f, %.2f. \n", */ 
	    /* rank, leftbuff[0], leftbuff[1], rightbuff[0], rightbuff[1]); */ 

    /* printf("From process %d: my data are %lf and %lf. \n", rank, localin[0], localin[1]); */ 
    /* MPI_Finalize(); */ 
    if(rank == 0){
	printf("Final data: \n"); 
	for(int i=0; i<num_values; i++){
	    printf("%.4f, ", alldata[i]); 
	}
	printf("\n");
    }
    MPI_Finalize(); 
   
}


int read_input(const char *file_name, double **values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "r"))) {
        perror("Couldn't open input file");
        return -1;
    }
    int num_values;
    if (EOF == fscanf(file, "%d", &num_values)) {
        perror("Couldn't read element count from input file");
        return -1;
    }
    if (NULL == (*values = malloc(num_values * sizeof(double)))) {
        perror("Couldn't allocate memory for input");
        return -1;
    }
    for (int i=0; i<num_values; i++) {
        if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
	    perror("Couldn't read elements from input file");
	    return -1;
	    }
    }
    if (0 != fclose(file)){
        perror("Warning: couldn't close input file");
    }
    return num_values;
}


int write_output(char *file_name, const double *output, int num_values) {
    FILE *file;
    if (NULL == (file = fopen(file_name, "w"))) {
	perror("Couldn't open output file");
	return -1;
    }
    for (int i = 0; i < num_values; i++) {
	if (0 > fprintf(file, "%.4f ", output[i])) {
	    perror("Couldn't write to output file");
	}
    }
    if (0 > fprintf(file, "\n")) {
    	perror("Couldn't write to output file");
    }
    if (0 != fclose(file)) {
    	perror("Warning: couldn't close output file");
    }
    return 0;
}
