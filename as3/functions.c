#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int cal_dim(int a){
    int result = 0; 
    while(a >>= 1){
	result++; 
    }
    return result; 
}

int read_input(const char* file_name, int** array){
    FILE* file;
    if((file = fopen(file_name, "r")) == NULL){
        perror("Couldn't open input file. \n");
        return -1;
    }

    int N;
    if(EOF == fscanf(file, "%d", &N)){
        perror("Couldn't read the number of elements. \n");
        return -1;
    }
    
    *array = malloc(N*sizeof(int)); 
    if(*array == NULL){
        perror("Couldn't allocate memory for array. \n");
        return -1;
    }
    
    for(int i=0; i<N; i++){
        if(EOF == fscanf(file, "%d ", &(*array)[i])){
	    perror("Couldn't read elements from the input file. \n");
	    return -1;
	}
    }
    
    if(0 != fclose(file)){
        perror("Warning: couldn't close input file. \n");
    }
    return N;
}

void vis(const int* array, int num){
    for(int i=0; i<num; i++){
	printf("%d ", array[i]); 
    }
    printf("\n"); 
}

int cmpfunc(const void* a, const void* b){
    return (*(int*)a - *(int*)b); 
}

/* int main(int argc, char** argv){ */
/*     int* array; */
/*     int num = read_input(argv[1], &array); */
/*     printf("Data read! \n"); */ 
/*     vis(array, num); */
/*     qsort(array, num, sizeof(int), cmpfunc); */ 
/*     vis(array, num); */ 
/* } */
