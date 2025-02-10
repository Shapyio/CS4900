// cs4900
// Project 05
// Shapiy Sagiev
// w1140sxs
// Due 24 March 2023 before 10pm
// System = Fry
// Compiler syntax = ./compile.sh proj05
// Job Control File = proj05.batch
// Additional File  = test05.dat    <-- sample input file for A
//                    test05_m.dat  <-- sample input file for r
// Results file     = proj05.pdf
//
// using row-cyclic distribution scheme
// global_row=rank + local_row * ncpu
// Similar to dealing cards
//
// If you know the rank then
// local_row = (global_row-rank)/ncpu
// ...
// If you do not know the rank then
// the rank that hold this row is the one that makes
// (global_row-rank)%ncpu == 0


#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>

//#define DEBUG

int main(int argc, char *argv[]);
void usage(void);
void get_n_m(int argc, char *argv[], int nAndm[2]);
void define_n_local(int *n_local, int ncpu, int n);
void initializeA(float **A, int n, int *n_local, int rank, int ncpu, int argc, char *argv[]);
void create_reference_sol(float **r, int *n_local, int n, int m, int rank, int ncpu, int argc, char *argv[]);
void parallel_matrix_multiply(float **A, float **r, float **y, int n, int m, int *n_local, int rank, int ncpu);
void pivot(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu, int active_row);
void swap_max_index(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu, int active_row);
void swapLocal(int global1, int global2,float **A, float **y, int ncpu, int n, int m);
void parallelSwap(int global1, int global2,float **A, float **y, int ncpu, int n, int m);
void gaussian_elimination(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu);
void update_rows(float **A,float **y,int n,int m,int *n_local,int rank,int ncpu,int active_row);
void back_substitution(float **A, float **y, float **x, int n, int m, int rank, int ncpu);
void rms_error(float **x, float **r, int *n_local, int n, int m, int rank);
// helper functions
void printA(float **A, int n, int rank, int cpu);
void printr(float **r, int n, int m, int rank, int cpu);
void printx(float **x, int n, int m, int rank, int cpu);
void printy(float **y, int n, int m, int rank, int cpu);
/**********************************************************/
int main(int argc, char *argv[]) {

  //make some names
  //n - global matrix size
  //n_local # of row on my rank

  // global matrix names from the instructions
  //          but define each locally
  // A        matrix             global (n x n) local (n_local x n)
  // m        # of rhs vectors
  //          (m=10 for undergrads)
  // r        reference matrix   global (n x m) local (n_local x m)
  // y        rhs matrix         global (n x m) local (n_local x m)
  // x        solution matrix    global (n x m) local (n_local x m)
  // difsq    (x-r)^2                           local (1 x m)
  // epsilon  matrix error vector               local (1 x m) rank 0 only


  // *******************************************//
  // step 1 obtain n (and m if grad student)
  //
  int nAndm[2], n, m;
  get_n_m(argc, argv, nAndm);
  n=nAndm[0];                 // global
  m=nAndm[1];                 // global
  printf("Completed step 1!\n");
  // *******************************************//
  // Standard MPI initialization calls
  //
  int rank, ncpu;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
  // Enable MPI error handling for MPI_COMM_WORLD
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  // *******************************************//
  // step 2 calculate n_local for each rank
  //
  int n_local[ncpu];
  define_n_local(n_local, ncpu, n);
  // if (rank == 0) {
  //   for (int i=0; i<ncpu; i++) {
  //     printf("n_local[%i]=%i\n", i, n_local[i]);
  //   }
  // }
  printf("Rank %i completed step 2!\n", rank);
  // *******************************************//
  // step 3 use n_local and m to create local memory
  // A, r, y, x, difsq, and if rank 0 epsilon
  //
  float **A = (float **)malloc(sizeof(float*) * n_local[rank]);
  float **r = (float **)malloc(sizeof(float*) * n_local[rank]);
  float **y = (float **)malloc(sizeof(float*) * n_local[rank]);
  float **x = (float **)malloc(sizeof(float*) * n_local[rank]);
  int row;
  for (row=0; row<n_local[rank]; row++){
    A[row] = (float*)malloc(sizeof(float) * n);
    r[row] = (float*)malloc(sizeof(float) * m);
    y[row] = (float*)malloc(sizeof(float) * m);
    x[row] = (float*)malloc(sizeof(float) * m);
  }
  printf("Rank %i completed step 3!\n", rank);

  // *******************************************//
  // step 4 read in the data for A
  //
  initializeA(A, n, n_local, rank, ncpu, argc, argv);
  // printA(A, n, rank, ncpu);
  printf("Rank %i completed step 4!\n", rank);
  // *******************************************//
  // step 5 create reference solution
  //
  create_reference_sol(r,n_local,n,m,rank,ncpu,argc,argv);
  // printr(r, n, m, rank, ncpu);
  printf("Rank %i completed step 5!\n", rank);
  // *******************************************//
  // step 6 Calculate y (RHS)
  //
  MPI_Barrier(MPI_COMM_WORLD);
  parallel_matrix_multiply(A,r,y,n,m,n_local,rank,ncpu);
  //printy(y, n, m, rank, ncpu);
  double wtime;
  wtime = MPI_Wtime();
  printf("Rank %i completed step 6!\n", rank);
  // *******************************************//
  // step 7 Gaussian Elimination
  //        Functions used:
  //        pivot
  //           find_max_index
  //           swap_rows
  //        update_rows
  //           broadcast_row
  //
  gaussian_elimination(A,y,n,m,n_local,rank,ncpu);
  //printA(A, n, rank, ncpu);
  printf("Rank %i completed step 7!\n", rank);
  // *******************************************//
  // step 8 Back Substitution
  //
  back_substitution(A, y, x, n, m, rank, ncpu);
  printf("Rank %i completed step 8!\n", rank);
  // *******************************************//
  // step 9 Calculate error and print results
  //
  //printA(A, n, rank, ncpu);
  //printr(r, n, m, rank, ncpu);
  //printx(x, n, m, rank, ncpu);
  rms_error(x, r, n_local, n, m, rank);
  printf("Rank %i completed step 9!\n", rank);
  wtime = MPI_Wtime() - wtime;
  if ( rank==0)
    printf ( "n=%d  Elapsed wallclock time is %g\n", n, wtime );

  for (row=0; row<n_local[rank]; row++){
    free(A[row]);
    free(r[row]);
    free(y[row]);
    free(x[row]);
  }
  free(A);
  free(r);
  free(y);
  free(x);

  printf("[%d] ending .....................\n",rank);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

}

void get_n_m(int argc, char *argv[], int nAndm[2]){

  // error values
  // 0 normal state
  // 1 flag DNE
  // 2 n < 1 matrix size invalid
  // 3 error while reading matrix file
  // 4 error while reading rhs file

  int c, n=0, m=10;
  FILE *file;

  nAndm[0]=n;
  nAndm[1]=m;

  //  optind=1;
  while(1){

    c = getopt(argc, argv, "r:d:f:");

    if( c == -1) {
      break; // nothing left to parse in the command line
    }

    switch(c) {
    case 'r':
      // printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg);
      nAndm[0]=atoi(optarg);
      break;  // break out of switch stmt
    case 'd':
      // printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg);
      nAndm[0]=atoi(optarg);
      break;  // break out of switch stmt
    case 'f':
      // printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg);
      file = fopen(optarg, "r");
      if (file == NULL) {
        fprintf(stderr, "Erorr: Unable to open file %s\n", optarg);
        exit(1);
      }
      break;  // break out of switch stmt
    default: fprintf(stderr,"This flag DNE: %c\n", optopt);
      exit(1);
    }
  }

  printf("get_n_m done ...... \n");

#ifdef DEBUG
  printf("FROM %s: n=%d m=%d\n",__func__,n,m);
#endif

}


void define_n_local(int *n_local, int ncpu, int n){

  // Calculates which processor has what part of the global elements
  int min_rows = (int)(n / ncpu); // All processors have at least this many rows
  int remainder = (n % ncpu); // Remainder of above calculation

  // Calculate the number of rows for each processor (rank)
  for (int i = 0; i < ncpu; i++) {
    n_local[i] = min_rows;
    if (i < remainder) {
      n_local[i]++;
    }
  }

#ifdef DEBUG
  printf("FROM %s: \n",__func__);
#endif

}


void initializeA(float **A, int n, int *n_local, int rank, int ncpu, int argc, char *argv[]){
  srand(rank + 1); // Seed randomizer
  optind=1; // Reset parsing of command line options

  while(1) {

    char c = getopt(argc, argv, "r:d:f:");

    if( c == -1) {
      break; // nothing left to parse in the command line
    }

    FILE *file; // Initialize the file
    switch(c) {
    case 'r':
      // Fill matrix A with random values
      for (int i=0;i<n_local[rank];i++) { // Loop thru n_local rows
        for (int j=0; j<n; j++) { // Loop thru columns
          A[i][j] = rand(); // Assign random number
        }
      }
      break;  // break out of switch stmt
    case 'd':
      // Fill matrix A with diagonal counting up
      for (int j = 0; j < n_local[rank]; j++) { // Loop thru the cpu rank's n_local rows
        int global_row = rank + j * ncpu; // global_row = rank + local_row * ncpu
        for (int k=0; k<n; k++) { // Loop thru row of matrix
          if (k == global_row) {
            A[j][k] = global_row + 1; // Diagonal elements count up from 1
          } else {
            A[j][k] = 0; // Non-diagonal elements are zero
          }
        }
      }
      break;  // break out of switch stmt
    case 'f':
      // Read matrix from file input
      file = fopen(optarg, "r"); // Assuming optarg contains the filename
      if (file == NULL) {
        fprintf(stderr, "Error: Unable to open file %s\n", optarg);
        exit(1);
      }

      int matrix_size;
      fscanf(file, "%d", &matrix_size); // Read matrix size from file
      if (matrix_size != n) {
        fprintf(stderr, "Error: Matrix size specified in file doesn't match input size\n");
        exit(1);
      }

      for (int i = 0; i < n_local[rank]; i++) {
        for (int j = 0; j < n; j++) {
          if (fscanf(file, "%f", &A[i][j]) != 1) {
            fprintf(stderr, "Error reading data from file\n");
            exit(1);
          }
        }
      }

      fclose(file); // Close file stream
      break;  // break out of switch stmt
    default: fprintf(stderr,"This flag DNE: %c\n", optopt);
      exit(1);
    }
  }
  // printf("initializeA done ...... \n");

#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}


void create_reference_sol(float **r, int *n_local, int n, int m, int rank, int ncpu, int argc, char *argv[]){

  for (int i=0; i<n_local[rank];i++) { // Loop thru n_local size of r (rows)
    for (int j=0; j<n; j++) { // Loop thru n size of r (columns)
      r[i][j] = rand(); // Set numbers as rand()
    }
  }

#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void MPI_Result_Check(int rank, int result, char *func) {
  if (result != MPI_SUCCESS) {
    // Handle error
    printf("***Rank %i, running %s, failed with error code: %d\n", rank, func, result);
  }
}

void parallel_matrix_multiply(float **A, float **r, float **y, int n, int m, int *n_local, int rank, int ncpu){

  // Trying to do
  //
  // y=A*x
  //
  // Watch video lecture from March 07
  // put your code here

  //    setting things up    //
  //  |                   |
  // \ /                 \ /

  // MPI variables
  MPI_Request request;
  int flag;
  MPI_Status status;

  // Allocate memory for compressed vector
  float *compressed_r = (float *)malloc(n_local[rank] * m * sizeof(float)); //FIXME: 

  // initialize y to zero
  for (int i = 0; i < n_local[rank]; i++) {
    for (int j = 0; j < m; j++) {
      y[i][j] = 0.0;
    }
  }
  // compress r matrix into a vector
  for (int i = 0; i < n_local[rank]; i++) {
    for (int j = 0; j < m; j++) {
      compressed_r[i * m + j] = r[i][j];
    }
  }

  //initialize
  int pass_to_rank=(rank+1) % ncpu;          // does not change in loop
  int get_from_rank=(rank-1+ncpu) % ncpu;         // does not change in loop
  int data_started_on_rank=rank;  // changes in loop

  // / \                 / \
  //  |                   |
  ///////////////////////////

  // Loop over all ranks and use this as a tag
  for (int count_cpu=0; count_cpu<ncpu; ++count_cpu){

    // Send compressed vector to neighboring rank (pass_to_rank)
    printf("\nLoop: %i, Rank %i sending %i to Rank %i, ", count_cpu, rank, n_local[data_started_on_rank] * m, pass_to_rank);
    for (int x=0; x<n_local[data_started_on_rank];x++) {
      printf("%f, ", compressed_r[x]);
    }
    int result = MPI_Isend(compressed_r, n_local[data_started_on_rank] * m, MPI_FLOAT, pass_to_rank, count_cpu, MPI_COMM_WORLD, &request);
    MPI_Result_Check(rank, result, "MPI_Isend");
    //unpack
    float **temp_r = (float **)malloc(sizeof(float*) * n_local[data_started_on_rank]); // Temp r
    for (int row=0; row<n_local[data_started_on_rank]; row++){
      temp_r[row] = (float*)malloc(sizeof(float) * m);
    }
    // Decompress compressed_r back into a matrix
    for (int i = 0; i < n_local[data_started_on_rank]; i++) {
      for (int j = 0; j < m; j++) {
        temp_r[i][j] = compressed_r[i * m + j];
      }
    }

    // multiply part from each rank filling up y matrix
    for (int i = 0; i < n_local[data_started_on_rank]; i++) {
      for (int j = 0; j < m; j++) {
        for (int k = 0; k < n; k++) {
          y[i][j] += A[i][k] * temp_r[k][j]; // FIXME: A = nxn temp_r = nx10
        }
      }
    }

    // update data_started_on_rank value
    --data_started_on_rank;
    if (data_started_on_rank < 0)
      data_started_on_rank=ncpu-1;

    MPI_Test(&request,&flag,&status);

    if (flag==1){ // MPI_Isend complete
      // Get compressed vector from get_from_rank
      result=MPI_Recv(compressed_r, n_local[data_started_on_rank] * m, MPI_FLOAT, get_from_rank, count_cpu, MPI_COMM_WORLD, &status);
      MPI_Result_Check(rank, result, "MPI_Recv");
    } else { //MPI_Isend not complete
      // Allocate memory for TEMPORARY compressed vector
      float *temp_compressed_r = (float *)malloc(n_local[data_started_on_rank] * m * sizeof(float));
      result=MPI_Recv(temp_compressed_r,m*n_local[data_started_on_rank],MPI_FLOAT,get_from_rank,count_cpu,MPI_COMM_WORLD,&status);
      MPI_Result_Check(rank, result, "MPI_Recv");
      MPI_Wait(&request, &status);        // this blocks until MPI_Isend is done
      for (int i=0; i<(m*n_local[data_started_on_rank]); ++i){
        compressed_r[i]=temp_compressed_r[i];
      }
      free(temp_compressed_r); // Free allocated memory
    }
    printf("\nLoop: %i, Rank %i receiving %i from Rank %i, ", count_cpu, rank, n_local[data_started_on_rank] * m, get_from_rank);
    for (int x=0; x<n_local[data_started_on_rank];x++) {
      printf("%f, ", compressed_r[x]);
    }
    for (int row=0; row<n_local[data_started_on_rank]; row++){
      free(temp_r[row]);
    }
    free(temp_r);
  }
  MPI_Barrier(MPI_COMM_WORLD); // Wait for all processors to finish
  // Free allocated memory
  free(compressed_r);

#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}


void gaussian_elimination(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu){

  // main loop can not be parallized
  int row;
  for (row=0; row<(n-1); row++){   // global loop - all ranks do this loop
    /*if (rank==0){
      printf("pivot starting row=%d of %d\n",row,(n-2));
    }
    printA(A, n, rank, ncpu); */
    pivot(A,y,n,m,n_local,rank,ncpu,row);
    /*     if (rank==0){
      printf("\n\n\nAfter pivot row=%d\n",row);
    }
    printA(A, n, rank, ncpu);
    printy(y, n, m, rank, ncpu);*/

    update_rows(A,y,n,m,n_local,rank,ncpu,row);
    /*if (rank==0){
      printf("\nAfter update_rows\n");

    }
    printA(A, n, rank, ncpu);
    printy(y, n, m, rank, ncpu); */

  }
  /*  if (rank==0){
    printf("pivot ending ..................\n");
  }

  printA(A, n, rank, ncpu);
  fflush(stdout); */

#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void pivot(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu, int active_row){

  swap_max_index(A,y,n,m,n_local,rank,ncpu,active_row);

#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void swap_max_index(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu, int active_row){
  printf("\n-RANK %i entering swap_max_index!", rank);
  struct {
    float number;
    int   index;
  } localMax, globalMax;
  printf("\n-RANK %i defined struct!", rank);
  localMax.number = A[active_row][active_row];
  localMax.index = active_row;
  printf("\n-RANK %i set initial localMax!", rank);
  for (int i = active_row + 1; i < n; i++) {
    if (A[i][active_row] > localMax.number) {
      localMax.number = A[i][active_row];
      localMax.index = i;
    }
  }
  printf("\n-RANK %i, got local max!", rank);
  MPI_Allreduce(&localMax, &globalMax, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  printf("\n-RANK %i All reduced!", rank);
  if (globalMax.index != active_row) {
    if (rank == globalMax.index % ncpu) {
      swapLocal(globalMax.index, active_row, A, y, ncpu, n, m);
    } else if (rank == active_row % ncpu) {
      parallelSwap(globalMax.index, active_row, A, y, ncpu, n, m);
    }
  }
  printf("\n-RANK %i, completed pivot math!", rank);
  // First find the max index and then use one of the below functions to swap if needed.
  // Note that the row being worked on "might" be the max value,... thus no swapping
  // Needed functions swapLocal and parallelSwap

#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void swapLocal(int global1, int global2, float **A, float **y, int ncpu, int n, int m){

  // Global row to local row math
  int local_row1 = (global1 - global1 % ncpu) / ncpu;
  int local_row2 = (global2 - global2 % ncpu) / ncpu;
  // Swap A
  for (int j = 0; j < n; j++) {
    float temp = A[local_row1][j];
    A[local_row1][j] = A[local_row2][j];
    A[local_row2][j] = temp;
  }

  // Swap y
  for (int j = 0; j < m; j++) {
    float temp = y[local_row1][j];
    y[local_row1][j] = y[local_row2][j];
    y[local_row2][j] = temp;
  }

#ifdef DEBUG
  printf("[%d]FROM %s: \n",global1%ncpu,__func__);
#endif

}


void parallelSwap(int global1, int global2,float **A, float **y, int ncpu, int n, int m){
  // Broadcast global1 and global2 to all ranks
  MPI_Bcast(&global1, 1, MPI_INT, global1 % ncpu, MPI_COMM_WORLD);
  MPI_Bcast(&global2, 1, MPI_INT, global2 % ncpu, MPI_COMM_WORLD);

  // Global to local row math
  int local_row1 = (global1 - global1 % ncpu) / ncpu;
  int local_row2 = (global2 - global2 % ncpu) / ncpu;
  // Swap A and Y
  for (int j = 0; j < n; j++) {
    float temp = A[local_row1][j];
      A[local_row1][j] = A[local_row2][j];
      A[local_row2][j] = temp;
  }

  for (int j = 0; j < m; j++) {
    float temp = y[local_row1][j];
    y[local_row1][j] = y[local_row2][j];
    y[local_row2][j] = temp;
  }

#ifdef DEBUG
  printf("[%d]FROM %s: \n",global1%ncpu,__func__);
#endif

}


void update_rows(float **A,float **y,int n,int m,int *n_local,int rank,int ncpu,int active_row){

  int local_n = *n_local;
  for (int i = 0; i < local_n; i++) {
    if (i != active_row) {
      float factor = A[i][active_row] / A[active_row][active_row];
      for (int j = active_row; j < m; j++) {
        A[i][j] -= factor * A[active_row][j];
      }
      y[i][0] -= factor * y[active_row][0];
    }
  }


#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void back_substitution(float **A, float **y, float **x, int n, int m, int rank, int ncpu){
  /*for  i <-- n-1 down to 1 do
     x[i] <-- b[i] / a[i,i]
     for j <-- 0 to i-1 do
         b[j] <-- b[j] - x[i] x a[j,i]
         a[j,i] <-- 0
     end for loop j
   end for loop i
  */

  for (int i=n-1; i>0; i--) {
    (*x)[i] = (*y)[i] / A[i][i];
    for (int j=0; j<i-1; j++) {
      (*y)[j] = (*y)[j] - (*x)[i] * A[j][i];
      A[j][i] = 0;
    }
  }

#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void rms_error(float **x, float **r, int *n_local, int n, int m, int rank){
  float epsilon[m], difsq[m];
  int rhs, row;
  for(rhs=0; rhs<m; ++rhs){
    difsq[rhs]=0.0;
    for(row=0; row<n_local[rank]; ++row)
      difsq[rhs]=difsq[rhs]+(x[row][rhs]-r[row][rhs])*(x[row][rhs]-r[row][rhs]);
  }
  MPI_Reduce(difsq, epsilon, m, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD );

  if (rank==0){
    for(rhs=0; rhs<m; rhs++){
      epsilon[rhs]=sqrt(epsilon[rhs])/n;
      printf("rhs=%02d epsilon=%f\n",rhs,epsilon[rhs]);
    }
  }
#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void printA(float **A, int n, int rank, int ncpu){
  int row, row_local,col;

  if(rank<1)
    printf("Matrix A dimensions %dx%d\n",n,n);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  for( row=0; row<n; ++row){
    row_local=row/ncpu;
    if ((row%ncpu)==rank){
      printf("[%d] ",rank);
      for( col=0; col<n; ++col){
        printf("%4.2f ",A[row_local][col]);
      }
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


void printr(float **r, int n, int m, int rank, int ncpu){
  int row, row_local,col;

  if(rank<1)
    printf("Matrix r dimensions %dx%d\n",n,m);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  for( row=0; row<n; ++row){
    row_local=row/ncpu;
    if ((row%ncpu)==rank){
      printf("[%d] ",rank);
      for( col=0; col<m; ++col){
        printf("%4.2f ",r[row_local][col]);
      }
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


void printx(float **x, int n, int m, int rank, int ncpu){
  int row, row_local,col;

  if(rank<1)
    printf("Solution matrix x dimensions %dx%d\n",n,m);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  for( row=0; row<n; ++row){
    row_local=row/ncpu;
    if ((row%ncpu)==rank){
      printf("[%d] ",rank);
      for( col=0; col<m; ++col){
        printf("%f ",x[row_local][col]);
      }
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void printy(float **y, int n, int m, int rank, int ncpu){
  int row, row_local,col;

  if(rank<1)
    printf("RHS matrix y dimensions %dx%d\n",n,m);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  for( row=0; row<n; ++row){
    row_local=row/ncpu;
    if ((row%ncpu)==rank){
      printf("[%d] ",rank);
      for( col=0; col<m; ++col){
        printf("%f ",y[row_local][col]);
      }
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
