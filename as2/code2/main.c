#include "mat.h"
#include <math.h>
#include <mpi.h>
#include <string.h>

#define A_TAG 1349
#define B_TAG 204
#define T_TAG 111

#define VIS 0
#define SAVE 0
#define TIMEANALYSIS 1

int main(int argc, char** argv){
    if(argc != 3){
	printf("Error: invalid parameters. Input file name and output file name required. \n"); 
	return -1; 
    }
    /* MPI data declaration */
    int rank, size; 
    MPI_Comm GRID_COMM; 
    MPI_Datatype checkboard_type;

    int grid_dim[2], grid_period[2], my_coord[2] = {0,0}, 
	up_coord[2], down_coord[2], 
	left, right, up, down; 
    char* input_file = argv[1]; 
    char* output_file = argv[2];

    /* MPI initialization */
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    MPI_Request Arequest[size], Brequest[size],Trequest, Trequestb, Trequestm;
    MPI_Status status[size]; 

    grid_dim[0] = grid_dim[1] = (int)sqrt((double)size); 
    grid_period[0] = grid_period[1] = 1; 
    
    /* if(rank == 0){printf("The grid dim are %d, %d. \n", grid_dim[0], grid_dim[1]); } */
    if(grid_dim[0]*grid_dim[1] != size){
	printf("Error: the number of cores cannot be grided. \n"); 
	MPI_Abort(MPI_COMM_WORLD, 1); 
	return -1; 
    }
    
    /* Create topology */
    MPI_Cart_create(MPI_COMM_WORLD, 2, grid_dim, grid_period, 1, &GRID_COMM); 
    MPI_Cart_coords(GRID_COMM, rank, grid_dim[0], my_coord);
    
    /* Find neighbours */

    up_coord[1] = down_coord[1] = my_coord[1]; 
    up_coord[0] = (my_coord[0] - 1 + grid_dim[0]) % grid_dim[0]; 
    down_coord[0] = (my_coord[0] + 1) % grid_dim[0]; 

    MPI_Cart_rank(GRID_COMM, up_coord, &up); 
    MPI_Cart_rank(GRID_COMM, down_coord, &down);

    /* Local variable declaration */
    int dim_all, dim; 
    float *A_all, *B_all, *C_all, *A_local, *A_cal, *B_cal, *C; 
    if(rank == 0){
	dim_all = read_input(argv[1], &A_all, &B_all); 
	if(dim_all <= 0){
	    printf("Error: Read input data failed! \n"); 
	    return -1; 
	}
	C_all = malloc(sizeof(float)*dim_all*dim_all); 
	/* printf("Original data dim is: %d. \n", dim); */
	dim = dim_all/(int)sqrt(size); 
	if(dim * (int)sqrt(size) != dim_all){
	    printf("The matrix cannot be grided. \n"); 
	    MPI_Abort(GRID_COMM, 1); 
	    return -1; 
	}
    }
    
    MPI_Bcast(&dim, 1, MPI_INT, 0, GRID_COMM); 
    /* printf("From process %d: local data dim is: %d. \n", rank, dim); */   
    
    /* Define checkboard data type (for task assigning) */
    /* MPI_Type_vector(dim, dim,dim*grid_dim[0], MPI_FLOAT, &checkboard_type); */ 
    /* MPI_Type_commit(&checkboard_type); */ 
    
    /* Another approach - cyclic checkboard split. */
    int block_len[dim]; 
    MPI_Aint indices[dim]; 
    MPI_Datatype data_type[dim];
     
    for(int i=0; i<dim; i++){
	block_len[i] = 1; 
	indices[i] = i*grid_dim[0]*sizeof(float);
	data_type[i] = MPI_FLOAT; 
    }
    /* indices[dim-1] = grid_dim[0]*dim*grid_dim[0]*sizeof(float); */ 
    /* if(rank == 0){ */
	/* printf("Dim is %d. \n===========================\n", dim); */ 
	/* for(int i=0; i<dim; i++){ */
	    /* printf("Index is %d. \n", indices[i]); */ 
	/* } */
    /* } */
    MPI_Type_create_struct(dim, block_len, indices, data_type, &checkboard_type); 
    MPI_Type_create_resized(checkboard_type, 0, grid_dim[0]*grid_dim[0]*dim*sizeof(float), &checkboard_type);
    MPI_Type_commit(&checkboard_type); 
    

    /* Allocate memory for local stored partial-A and calculation-required partial A and B */
    A_local = malloc(dim*dim*sizeof(float)); 
    A_cal = malloc(dim*dim*sizeof(float)); 
    B_cal = malloc(dim*dim*sizeof(float)); 
    C = calloc(dim*dim, sizeof(float)); 

    if(A_local == NULL || A_cal == NULL || B_cal == NULL || C == NULL){
	printf("From rank %d: memory allocation for matrices failed! \n", rank); 
	return -1; 
    }


    /* Assign matrix A and B to all processes */ 
    if(rank == 0){
	for(int i=0; i<size; i++){
	    /* printf("Index is: %d. \n", i%grid_dim[0]*dim*grid_dim[0] + i/grid_dim[0]); */ 
	    MPI_Isend(&(A_all[(i%grid_dim[0]) + i/grid_dim[0]*dim*grid_dim[0]]),  
		    dim, checkboard_type, i, A_TAG+i, GRID_COMM, &Arequest[i]); 
	    MPI_Isend(&(B_all[(i%grid_dim[0]) + i/grid_dim[0]*dim*grid_dim[0]]),  
		    dim, checkboard_type, i, B_TAG+i, GRID_COMM, &Brequest[i]); 
	}
    }
    MPI_Recv(A_local, dim*dim, MPI_FLOAT, 0, A_TAG+rank, GRID_COMM, &status[0]); 
    MPI_Recv(B_cal, dim*dim, MPI_FLOAT, 0, B_TAG+rank, GRID_COMM, &status[0]); 
    if(rank == 0){
	MPI_Waitall(size, Arequest, status); 
	MPI_Waitall(size, Brequest, status); 
    }
    /* if(rank == 0){ */
	/* vis(A_local, dim); */ 
    /* } */ 
    /* Define Row communicator */
    MPI_Comm ROW_COMM;  
    MPI_Comm_split(GRID_COMM, my_coord[0], rank, &ROW_COMM); 
#if TIMEANALYSIS
	/* Communication time included */
	double start, my_execution_time,m_start,m_cost, b_start, b_cost;
	double* exe_time;
	double* b_time;
	double* m_time;
	//double* ser_time;
#endif
    /* Calculation for A*B parallel*/
    for(int i=0; i<grid_dim[0]; i++){
#if TIMEANALYSIS
		start = MPI_Wtime();
#endif
	/* Broadcast target partial matrix A */
	m_start = MPI_Wtime();
	if((my_coord[0]+i)%grid_dim[0] == my_coord[1]){
	    memcpy(A_cal, A_local, sizeof(float)*dim*dim); 
	}
	MPI_Bcast(A_cal, dim*dim, MPI_FLOAT, 
		(my_coord[0]+i)%grid_dim[0], ROW_COMM);
	m_cost += MPI_Wtime() - m_start;
	/* Local matrix calculation */
	/* printf("Process %d: round %d: My A is %.3f, my B is %.3f. \n", rank, i, A_cal[0], B_cal[0]); */ 

	//ser_start = MPI_Wtime();
	float* temp = multiply(A_cal, B_cal, dim); 
	//ser_cost += MPI_Wtime()-start;
	sum(C, temp, dim);
	free(temp);
	m_start = MPI_Wtime();
	MPI_Isend(B_cal, dim*dim, MPI_FLOAT, up, B_TAG+rank, GRID_COMM, &Brequest[0]); 
	MPI_Recv(B_cal, dim*dim, MPI_FLOAT, down, B_TAG+down, GRID_COMM, &status[0]); 
	MPI_Wait(&Brequest[0], &status[0]);
	m_cost += MPI_Wtime() - m_start;
	b_start = MPI_Wtime();
	MPI_Barrier(GRID_COMM);
	b_cost += MPI_Wtime() - b_start;
#if TIMEANALYSIS
	/* Time evaluation part */
	my_execution_time += MPI_Wtime() - start;
#endif
    }

    /* Collect result from all processes to process with rank 0 */ 
    MPI_Isend(C, dim*dim, MPI_FLOAT, 0, B_TAG+rank, GRID_COMM, &Brequest[0]); 
#if TIMEANALYSIS
	MPI_Isend(&my_execution_time, 1, MPI_DOUBLE, 0, T_TAG* rank, GRID_COMM, &Trequest);
	MPI_Isend(&b_cost, 1, MPI_DOUBLE, 0, 7777* rank, GRID_COMM, &Trequestb);
	MPI_Isend(&m_cost, 1, MPI_DOUBLE, 0, 6666* rank, GRID_COMM, &Trequestm);
	//MPI_Isend(&ser_cost, 1, MPI_DOUBLE, 0, 777* rank, GRID_COMM, &Trequests);
#endif
    if(rank == 0){
#if TIMEANALYSIS
		exe_time = malloc(size * sizeof(double));
		b_time = malloc(size * sizeof(double));
		m_time = malloc(size * sizeof(double));
		//ser_time = malloc(size * sizeof(double));
		exe_time[0] = my_execution_time;
		b_time[0] = b_cost;
		m_time[0] = m_cost;
		//ser_time[0] = ser_cost;
#endif
	for(int i=0; i<size; i++){
	    MPI_Recv(&(C_all[(i%grid_dim[0]) + i/grid_dim[0]*dim*grid_dim[0]]), 
		    dim, checkboard_type, i, B_TAG+i, GRID_COMM, &status[0]); 
#if TIMEANALYSIS
		MPI_Recv(&exe_time[i], 1, MPI_DOUBLE, i, T_TAG* i, GRID_COMM, &status[0]);
		MPI_Recv(&b_time[i], 1, MPI_DOUBLE, i, 7777* i, GRID_COMM, &status[0]);
		MPI_Recv(&m_time[i], 1, MPI_DOUBLE, i, 6666* i, GRID_COMM, &status[0]);
		//MPI_Recv(&ser_time[i], 1, MPI_DOUBLE, i, 777* i, GRID_COMM, &status[0]);
#endif
	}

#if TIMEANALYSIS
	double max_time = 0, max_serial_time = 0,mm=0,mb=0;
	for (int i = 0; i < size; i++) {
		if (exe_time[i] > max_time) {
			max_time = exe_time[i];
		}
		if (b_time[i] > mb) {
			mb = b_time[i];
		}
		if (m_time[i] >mm) {
			mm = m_time[i];
		}
		/*if (ser_time[i] > max_serial_time) {
			max_serial_time = ser_time[i];
		}*/
	}
	printf("%f\n", max_time);
	printf("%f\n", mb);
	printf("%f\n", mm);
	//printf("%f\n", max_serial_time);
#endif

#if VIS
	/* Print result out */
	printf("Matrix A: \n"); 
	vis(A_all, dim*grid_dim[0]); 
	printf("Matrix B: \n"); 
	vis(B_all, dim*grid_dim[0]); 
	printf("A * B: \n"); 
	vis(C_all, dim*grid_dim[0]);
#endif
#if SAVE
	if(write_output(output_file, C_all, dim*grid_dim[0]) != 0){
	    printf("Write output matrix failed! \n"); 
	}
#endif
#if TIMEANALYSIS
	MPI_Wait(&Trequest, &status[0]);
	MPI_Wait(&Trequestb, &status[0]);
	MPI_Wait(&Trequestm, &status[0]);
	//MPI_Wait(&Trequests, &status[0]);
	free(exe_time);
	//free(ser_time);
#endif
	free(C_all); 
	free(A_all); 
	free(B_all); 
    }
    MPI_Wait(&Brequest[0], &status[0]); 
    free(A_local); 
    free(A_cal); 
    free(B_cal); 
    free(C); 
    MPI_Type_free(&checkboard_type); 

    MPI_Finalize();
    return 0;  
}

/*
const double P[15*7]= {
	1,0,0,0,0,0,0,
	-1,0,0,0,0,0,0,
	-1,0,1,0,0,0,0,
	0,1,0,0,0,0,0,
	0,-1,0,0,0,0,0,
	0,-1,0,1,0,0,0,
	0,0,-1,0,0,0,0,
	0,0,0,-1,0,1,0,
	0,0,0,0,-1,0,0,
	0,0,0,0,-1,0,1,
	0,0,0,0,0,-1,0,
	1,0,0,0,0,-1,0,
	0,0,0,0,0,0,-1};
*/