#include <stdio.h>
#include <stdlib.h>
#include <time.h>
int write_output(char*, const int*, int); 

int main(int argc, char** argv){
    int dim = atoi(argv[1]); 
    int matrix[dim*dim*2]; 
    for(int i=0; i<dim*dim*2; i++){
	/* srand(i); */ 
	matrix[i] = rand()%20;
	printf("%d \n", matrix[i]); 
    }
    char output[] = "data"; 
    write_output(output, matrix, dim); 
    return 0; 
}

int write_output(char* file_name, const int* output, int dim){
    FILE* file; 
    int num_values = dim*dim*2; 
    if (NULL == (file = fopen(file_name, "w"))) {
	    perror("Couldn't open output file");
		return -1;
    }

    fprintf(file, "%d ", dim); 

    for(int i=0; i<num_values; i++){
	if(0>fprintf(file, "%d ", output[i])){
	    perror("Couldn't write to output file");
	}
    }

    if(0>fprintf(file, "\n")){
	perror("Couldn't write to output file");
    }

    if(0!=fclose(file)){
        perror("Warning: couldn't close output file");
    }
    return 0;
}
