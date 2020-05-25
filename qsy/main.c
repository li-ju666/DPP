#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define COLLECTTAG 111
#define TIMETAG 204

#define TIMEANALYSIS 1
#define SAVE 0
#define HISTOGRAM 1

void sum_a0(double *w, double R, double *a0);
int find_r(double *w, double a0u2, double R);
void prop(int *x, double *w);
int write_output(char* file_name, int* output, int num_values);

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Usage: number_of_experiments output_file \n");
		return 1;
	}
	MPI_Status status;
	MPI_Request srequest, rrequest;
	
	int num_steps = atoi(argv[1]);
	char* output_name = argv[2];
	int rank, size, num_experiment;
	int *experiment_result;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		num_experiment = num_steps / size;

		if (num_experiment * size != num_steps) {
			printf("Error: the number of experiment cannot be divided. \n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return -1;
		}
	}
	MPI_Bcast(&num_experiment, 1, MPI_INT, 0, MPI_COMM_WORLD);
	const int x0[7] = { 900, 900, 30, 330, 50, 270, 20 };
	const double T = 100;
	const int m = 15;
	const int n = 7;
	const int R = 15;
	const int b = 20;
	int * P;
	P = (int *)calloc(m * n, sizeof(int *));
	P[0 * 7 + 0] = 1;
	P[1 * 7 + 0] = -1;
	P[2 * 7 + 0] = -1;
	P[2 * 7 + 2] = 1;
	P[3 * 7 + 1] = 1;
	P[4 * 7 + 1] = -1;
	P[5 * 7 + 1] = -1;
	P[5 * 7 + 3] = 1;
	P[6 * 7 + 2] = -1;
	P[7 * 7 + 2] = -1;
	P[7 * 7 + 4] = 1;
	P[8 * 7 + 3] = -1;
	P[9 * 7 + 3] = -1;
	P[9 * 7 + 5] = 1;
	P[10 * 7 + 4] = -1;
	P[11 * 7 + 4] = -1;
	P[11 * 7 + 6] = 1;
	P[12 * 7 + 5] = -1;
	P[13 * 7 + 0] = 1;
	P[13 * 7 + 6] = -1;
	P[14 * 7 + 6] = -1;
	int r;
	int *x, *result;
	double *w;
	double t, u1, u2, a0, delta_t, a0u2;
	x = (int *)malloc(n * sizeof(int));
	w = (double *)malloc(m * sizeof(double));
	result = (int *)malloc(n * num_experiment * sizeof(int));
#if TIMEANALYSIS
	double start, my_execution_time, exe_time;
	start = MPI_Wtime();
#endif
	for (int col = 0; col < num_experiment; col++) {
		t = 0;
		for (int i = 0; i < n; i++) {
			x[i] = x0[i];
		}
		while (t < T) {
			a0 = 0;
			prop(x, w);
			sum_a0(w, R, &a0);
			//printf("a0=%lf\n",a0);
			u1 = ((double)rand()) / RAND_MAX;
			u2 = ((double)rand()) / RAND_MAX;
			//printf("u1=%lf, u2 = %lf\n",u1,u2);
			delta_t = -log(u1) / a0;
			//printf("deltat=%lf\n",delta_t);
			a0u2 = u2 * a0;
			r = find_r(w, a0u2, R);
			//	printf("r=%d\n",r);
			for (int i = 0; i < n; i++) {
				x[i] += P[r * n + i];
			}
			t += delta_t;
		}
		for (int i = 0; i < n; i++) {
			result[col * n + i] = x[i];
		}
	}

#if HISTOGRAM
	int max_num = 0, min_num = 10000, total_max, total_min;
	double interval;
	double *interval_value;
	int * interval_num, *total_num;
	interval_value = (double *)malloc((b + 1) * sizeof(double));
	interval_num = (int *)calloc(b, sizeof(int));
	for (int i = 0; i < num_experiment; i++) {
		if (result[i * n] > max_num) {
			max_num = result[i * n];
		}
		if (result[i * n] < min_num) {
			min_num = result[i * n];
		}
	}
	MPI_Reduce(&max_num, &total_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&min_num, &total_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	
	if (rank == 0) {
		interval = (total_max - total_min) / b;
		for (int i = 0; i <=b; i++) {
			interval_value[i] = total_min + i * interval;
			printf("%lf, ", interval_value[i]);	
		}
		printf("\n");
	}
	MPI_Bcast(interval_value, (b+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (int i = 0; i < b; i++) {
		for (int j = 0; j < num_experiment; j++) {
			if (result[j*n] <= interval_value[i + 1] && result[j*n] > interval_value[i]) {
				interval_num[i] += 1;
			}
		}
	}
	//if (rank == 0) {
	total_num = malloc(b * sizeof(int));
	
	for (int i = 0; i < b; i++) {
		MPI_Reduce(&(interval_num[i]),&(total_num[i]), 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if(rank == 0){		
			printf("%d, ", total_num[i]);}
	}
//if (rank == 0) {
	free(total_num);
	free(interval_num);
	free(interval_value);

#endif


#if TIMEANALYSIS
	my_execution_time = MPI_Wtime() - start;
#endif
	if (rank == 0) {
		experiment_result = (int *)malloc(num_steps * n * sizeof(int));
		memcpy(experiment_result, result, sizeof(int) * num_experiment * n);
		for (int i = 1; i < size; i++) {
			MPI_Recv(&(experiment_result[i * num_experiment * n]), num_experiment * n,
				MPI_INT, i, COLLECTTAG+i, MPI_COMM_WORLD, &status);//warning type cart needed

		}
	}
	else {
		MPI_Send(result, num_experiment * n, MPI_INT, 0, COLLECTTAG + rank, MPI_COMM_WORLD);
	}
#if TIMEANALYSIS	
	MPI_Reduce(&my_execution_time,&exe_time,1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		printf("\n%lf\n", exe_time);
	}
#endif
	free(x);
	free(w);
#if SAVE
	int *susceptible;
	if (rank == 0) {
		susceptible = (int *)malloc(num_steps * sizeof(int));
		for (int i = 0; i < num_steps; i++) {
			susceptible[i] = experiment_result[i*n];
		}
		write_output(output_name, susceptible, num_steps);//warning type cart needed
	}
	free(susceptible);
#endif
	free(result);
	MPI_Finalize();
}
