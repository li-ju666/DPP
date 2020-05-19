#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define FINAL_TIME 100
#define SEND_TAG 1349 

void prop(int*, double*);  
int find_r(double*, double); 
void vis(int*, int, int); 
int write_output(char*, const int*, int); 
    
int main(int argc, char** argv){
    if(argc != 3){
	printf("Error: Number of experiments on each processor required! \n"); 
	return -1; 
    }

    int size, rank; 
    MPI_Status status; 

    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

    int sim_num = atoi(argv[1]); 
    double T = FINAL_TIME; 
    double t = 0; 
    //int x[7] = {900, 900, 30, 330, 50, 270, 20}; 
    double w[15]; 
    double a0, u1, u2; 
    double tau; 
    int r; 
    int* P = calloc(sizeof(int), 7*15); 
    int* result = calloc(sizeof(int), sim_num*7); 
    int x0[7] = {900, 900, 30, 330, 50, 270, 20}; 
    int x[7];
    int* cases = malloc(sizeof(int)*sim_num);  
    P[0*7+0] = 1; 
    P[1*7+0] = -1; 
    P[2*7+0] = -1; 
    P[2*7+2] = 1; 
    P[3*7+1] = 1; 
    P[4*7+1] = -1; 
    P[5*7+1] = -1; 
    P[5*7+3] = 1; 
    P[6*7+2] = -1; 
    P[7*7+2] = -1; 
    P[7*7+4] = 1; 
    P[8*7+3] = -1; 
    P[9*7+3] = -1; 
    P[9*7+5] = 1; 
    P[10*7+4] = -1; 
    P[11*7+4] = -1; 
    P[11*7+6] = 1; 
    P[12*7+5] = -1; 
    P[13*7+0] = 1; 
    P[13*7+6] = -1; 
    P[14*7+6] = -1; 

    for(int i=0; i<sim_num; i++){
	/* printf("Now i is: %d. \n", i); */ 
	memcpy(x, x0, sizeof(int)*7); 
	
	/* Assign array x to be x0 state */
	t = 0; 	
	int count = 0; 	
	while(t<T){
	    /* Compute w = prob(x) */
	    prop(x, w); 
	    
	    /* Sum w to get a0 */
	    a0 = 0; 
	    for(int j=0; j<15; j++){
		a0 += w[j]; 
	    }
	    
	    /* generate random numbers u1 and u2 */
	    u1 = (float)rand()/RAND_MAX; 
	    u2 = (float)rand()/RAND_MAX;
	    
	    /* Calculate tau */
	    tau = -log(u1)/a0; 
	    r = find_r(w, a0*u2); 
	    
	    for(int j=0; j<7; j++){
		x[j] += P[r*7+j]; 
	    }
	    t = t + tau; 
	    count++; 
	}
	/* vis(x, 1, 7); */ 
	memcpy(&(result[i*7]), x, 7*sizeof(int)); 
    }
    //printf("Calculation finished! \n");
    //vis(result, sim_num, 7);  
    
    for(int i=0; i<sim_num; i++){
	cases[i] = result[i*7]; 
    }
    
    int* cases_all; 
    if(rank == 0){
	cases_all = malloc(sizeof(int)*sim_num*size); 
	memcpy(cases_all, cases, sizeof(int)*sim_num); 
	
	for(int i=1; i<size; i++){
	    MPI_Recv(&(cases_all[i*sim_num]), sim_num, MPI_INT, i, SEND_TAG+i, MPI_COMM_WORLD, &status); 
	}
	
	if(write_output(argv[2], cases, sim_num)<0){
	    printf("Failed to write result! \n"); 
	}
	/* vis(cases_all, 1, sim_num*size); */ 
	free(cases_all); 
    
    }
    else{
	/* printf("???"); */ 
	MPI_Send(cases, sim_num, MPI_INT, 0, SEND_TAG+rank, MPI_COMM_WORLD); 
    }

    free(P); 
    free(result);
    free(cases); 
    MPI_Finalize();  
}


int find_r(double* array, double value){
    double sum = 0; 
    for(int i=0; i<15; i++){
	sum += array[i]; 
	if(sum >= value){
	    /* printf("Index is: %d. \n", i); */ 
	    return i; 
	}
    }
}

void vis(int* mat, int row, int col){
    for(int i=0; i<row; i++){
	for(int j = 0; j<col; j++){
	    printf("%d ", mat[i*col+j]); 
	}
	printf("\n"); 
    }
}

int write_output(char* file_name, const int* output, int num){
    FILE* file; 
    int num_values = num; 
    if (NULL == (file = fopen(file_name, "w"))) {
	    perror("Couldn't open output file");
		return -1;
    }

    for(int i=0; i<num_values; i++){
	if(0>fprintf(file, "%d, ", output[i])){
	    perror("Couldn't write to output file");
	    return -1; 
	}
    }

    if(0>fprintf(file, "\n")){
	perror("Couldn't write to output file");
	return -1; 
    }
    if(0!=fclose(file)){
        perror("Warning: couldn't close output file");
	return -1; 
    }
    return 0;
}

/* int main(){ */
/*     double a[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}; */ 
/*     double value = 0.3; */ 
/*     return find_r(a, value); */ 
/* } */
