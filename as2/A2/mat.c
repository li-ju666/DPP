#include "mat.h"

/* Function for matrix multiplication */
float* multiply(float* A, float* B, int dim) {
	float* result = calloc(dim*dim, sizeof(float));
	float* row = calloc(dim, sizeof(float));
	for (int j = 0; j < dim; j++) {
		for (int index = 0; index < dim; index++) {
			row[index] = B[index*dim + j];
		}
		for (int i = 0; i < dim; i++) {
			for (int index = 0; index < dim; index++) {
				result[i*dim + j] += A[i*dim + index] * row[index];
			}
		}
	}
	free(row);
	return result;
}

void sum(float* A, float* B, int dim){
    for(int i=0; i<dim*dim; i++){
	A[i] += B[i]; 
    }
}

/* Function to print matrix */
void vis(float* A, int dim){
    for(int i=0; i<dim; i++){
	for(int j=0; j<dim; j++){
	    printf("%.2f ", A[i*dim+j]); 
	}
	printf("\n"); 
    }
//    printf("-------------------------------------------------\n"); 
}

/* Function to read data from input file */
int read_input(const char* file_name, float** A, float** B){
    FILE *file;
    if ((file = fopen(file_name, "r")) == NULL) {
        perror("Couldn't open input file");
        return -1;
    }

    int N;
    if (EOF == fscanf(file, "%d", &N)) {
        perror("Couldn't read dimension of matrices from input file");
        return -1;
    }

    if (NULL == (*A = malloc(N*N * sizeof(float)))) {
        perror("Couldn't allocate memory for Matrix A");
        return -1;
    }
    for (int i=0; i<N*N; i++) {
        if (EOF == fscanf(file, "%f", &((*A)[i]))) {
	    perror("Couldn't read elements for matrix A from input file");
	    return -1;
	    }
    }

    if (NULL == (*B = malloc(N*N * sizeof(float)))) {
        perror("Couldn't allocate memory for Matrix B");
        return -1;
    }
    for (int i=0; i<N*N; i++) {
        if (EOF == fscanf(file, "%f", &((*B)[i]))) {
	    perror("Couldn't read elements for matrix B from input file");
	    return -1;
	    }
    }

    if (0 != fclose(file)){
        perror("Warning: couldn't close input file");
    }
    return N;

}
int write_output(char* file_name, const float* output, int dim){
    FILE* file; 
    int num_values = dim*dim; 
    if (NULL == (file = fopen(file_name, "w"))) {
	    perror("Couldn't open output file");
		return -1;
    }

    for(int i=0; i<num_values; i++){
	if(0>fprintf(file, "%.6f ", output[i])){
	    perror("Couldn't write to output file");
	    return -1; 
	}
	if((i+1)%dim == 0 && i+1 < num_values){
	    fprintf(file, "\n"); 
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
/* int main(int argc, char** argv){ */
/*     float *A, *B; */ 
/*     int dim; */ 
/*     dim = read_input(argv[1], &A, &B); */ 
/*     vis(A, dim); */ 
/*     vis(B, dim); */ 
/* } */
