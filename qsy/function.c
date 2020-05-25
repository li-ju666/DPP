#include <stdio.h>
#include <stdlib.h>

void sum_a0(double *w,double R,double *a0){
	for (int i = 0; i < R; i++) {
		*a0 += w[i];
	}
}

int find_r(double *w,double a0u2, double R) {
	int i = 0;
	double value = 0;
	for(i = 0;i < R; i++) {
		value += w[i];
		if(a0u2 <= value){return i;}
	}
	printf("error!");
	exit(-1);
}


int write_output(char* file_name, int* output, int num_values) {
	FILE* file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}

	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%d ", output[i])) {
			perror("Couldn't write to output file");
			return -1;
		}
	}

	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
		return -1;
	}
	return 0;
}
