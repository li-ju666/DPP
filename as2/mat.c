#include "mat.h"

/* Function for matrix multiplication */
float* multiply(float* A, float* B, int dim){
    float* result = calloc(dim*dim, sizeof(float)); 
    for(int i=0; i<dim; i++){
	for(int j=0; j<dim; j++){
	    for(int index=0; index<dim; index++){
		result[j*dim+i] += A[index*dim+i]*B[j*dim+index]; 
	    }
	}
    }
    return result; 
}

void sum(float* A, float* B, int dim){
    for(int i=0; i<dim*dim; i++){
	A[i] += B[i]; 
    }
    free(B); 
}

/* Function to print matrix */
void vis(float* A, int dim){
    for(int i=0; i<dim; i++){
	for(int j=0; j<dim; j++){
	    printf("%.5f ", A[j*dim+i]); 
	}
	printf("\n"); 
    }
    printf("-------------------------------------------------\n"); 
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

/* int main(int argc, char** argv){ */
/*     float *A, *B; */ 
/*     int dim; */ 
/*     dim = read_input(argv[1], &A, &B); */ 
/*     vis(A, dim); */ 
/*     vis(B, dim); */ 
/* } */
