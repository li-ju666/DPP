#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define FINAL_TIME 100
#define SEND_TAG 1349 
#define HIST_NUM 20

void prop(int*, double*);  
int find_r(double*, double); 
void vis(int*, int, int); 
int write_output(char*, const int*, const int*, int); 
int min(int*, int); 
int max(int*, int); 
void hist_stat(int*, int, int*, int, int**); 

int main(int argc, char** argv){
    /* check format of input parameters */
    if(argc != 3){
	printf("Error: Number of experiments on each processor and target file name required! \n"); 
	return -1; 
    }
    
    /* preparation for MPI initialization */
    int size, rank; 
    MPI_Status status; 

    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

    /* Local data preparation */
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
    int hist_range[HIST_NUM+1], max_value, min_value; 
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
    
    /* double start = MPI_Wtime(); */     
    /* Main loop for calculation */
    for(int i=0; i<sim_num; i++){
	/* Assign initial x0 state */
	memcpy(x, x0, sizeof(int)*7); 
	
	t = 0; 	
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
	    
	    /* Update x vector */
	    for(int j=0; j<7; j++){
		x[j] += P[r*7+j]; 
	    }
	    t = t + tau; 
	}
	/* Save result to result matrix */
	memcpy(&(result[i*7]), x, 7*sizeof(int)); 
    }
   
    /* Collect local susceptible people data to array cases */
    for(int i=0; i<sim_num; i++){
	cases[i] = result[i*7]; 
    }
    
    /* Check local maximum and minimum values */
    max_value = max(cases, sim_num); 
    min_value = min(cases, sim_num); 
    
    /* Check global maximum and minimum value */
    MPI_Allreduce(&max_value, &(hist_range[HIST_NUM+1]), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    MPI_Allreduce(&min_value, hist_range, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); 
   
    /* Split histogram range with intervals */
    float global_max = (float)hist_range[HIST_NUM+1]; 
    float global_min = (float)hist_range[0]; 
    float interval = (global_max - global_min)/HIST_NUM; 
    for(int i=0; i<HIST_NUM+1; i++){
	hist_range[i] = (int)global_min+(int)(i*interval); 
    }

    /* Do statistics on as-obtained cases and save the histogram result in the array named hist_result */
    int* hist_result; 
    hist_stat(cases, sim_num, hist_range, HIST_NUM, &hist_result); 

    /* Collect local statistics to process with rank 0 */
    int* hist_all; 
    if(rank == 0){
	hist_all = malloc(sizeof(int)*HIST_NUM); 
    }
    MPI_Reduce(hist_result, hist_all, HIST_NUM, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    

    /* double tcost = MPI_Wtime() - start; */ 
    /* double global_t; */ 
    /* MPI_Reduce(&tcost, &global_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); */ 
    /* Save global statistics to output file */
    if(rank == 0){
	/* printf("Time cost: %f. \n", global_t); */ 
	if(write_output(argv[2], hist_range, hist_all, HIST_NUM+1)<0){
	    printf("Failed to write result! \n"); 
	}
	free(hist_all); 
    }

    /* Free resources */
    free(P); 
    free(result);
    free(cases); 
    free(hist_result); 

    MPI_Finalize();  
}

/* Function to find r in the algorithm */
int find_r(double* array, double value){
    double sum = 0; 
    for(int i=0; i<15; i++){
	sum += array[i]; 
	if(sum >= value){
	    return i; 
	}
    }
}

/* Function to calculate the maximum integer of a given array */
int max(int* array, int num){
    int max = 0; 
    for(int i=0; i<num; i++){
	if(array[i]>max){
	    max = array[i]; 
	}
    }
    return max; 
}

/* Function to calculate the minimum integer of a given array */
int min(int* array, int num){
    int min=10000; 
    for(int i=0; i<num; i++){
	if(array[i]<min){
	    min = array[i]; 
	}
    }
    return min; 
}

/* Function to do statistics on given data and histogram range */
void hist_stat(int* array, int num, int* range, int ran_num, int** result){
    *result = calloc(sizeof(int), ran_num); 
    int count=0; 
    for(int i=0; i<num; i++){
	for(int j=0; j<ran_num; j++){
	    if(array[i] >= range[j] && array[i] < range[j+1]){
		(*result)[j]++; 
		count++; 
		break; 
	    }
	}
    }
    (*result)[ran_num-1] += num-count; 
}

/* Function to visualize matrix */
void vis(int* mat, int row, int col){
    for(int i=0; i<row; i++){
	for(int j = 0; j<col; j++){
	    printf("%d ", mat[i*col+j]); 
	}
	printf("\n"); 
    }
}

/* Function to write output to file */
int write_output(char* file_name, const int* output1, const int* output2, int num){
    FILE* file; 
    int num_values = num;
    if (NULL == (file = fopen(file_name, "w"))) {
	    perror("Couldn't open output file");
		return -1;
    }

    for(int i=0; i<num_values; i++){
	if(0>fprintf(file, "%d, ", output1[i])){
	    perror("Couldn't write to output file");
	    return -1; 
	}
    }
    if(0>fprintf(file, "\n")){
	perror("Couldn't write to output file");
	return -1; 
    }
    num_values--; 
    for(int i=0; i<num_values; i++){
	if(0>fprintf(file, "%d, ", output2[i])){
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
    return num_values;
}

