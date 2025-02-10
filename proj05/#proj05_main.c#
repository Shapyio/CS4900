// cs4900
// Project 05
// John Nehrbass
// w006jwn
// Due 24 March 2023 before 10pm
// System = Fry or Owens
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

  // *******************************************//
  // Standard MPI initialization calls
  //
  int rank, ncpu;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

  // *******************************************//
  // step 2 calculate n_local for each rank
  //
  int n_local[ncpu];
  define_n_local(n_local, ncpu, n);
  
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


  // *******************************************//
  // step 4 read in the data for A
  //
  initializeA(A, n, n_local, rank, ncpu, argc, argv);
  //printA(A, n, rank, ncpu);

  // *******************************************//
  // step 5 create reference solution
  //
  create_reference_sol(r,n_local,n,m,rank,ncpu,argc,argv);
  //printr(r, n, m, rank, ncpu);

  // *******************************************//
  // step 6 Calculate y (RHS)
  //
  parallel_matrix_multiply(A,r,y,n,m,n_local,rank,ncpu);
  //printy(y, n, m, rank, ncpu);
  double wtime;
  wtime = MPI_Wtime ( );

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

  // *******************************************//
  // step 8 Back Substitution
  //
  back_substitution(A, y, x, n, m, rank, ncpu);

  // *******************************************//
  // step 9 Calculate error and print results
  //
  //printA(A, n, rank, ncpu);
  //printr(r, n, m, rank, ncpu);
  //printx(x, n, m, rank, ncpu);
  rms_error(x, r, n_local, n, m, rank);

  wtime = MPI_Wtime ( ) - wtime;
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

  // put your code here
  
  int n,m;
  
  nAndm[0]=n;
  nAndm[1]=m;

#ifdef DEBUG
  printf("FROM %s: n=%d m=%d\n",__func__,n,m);
#endif

}


void define_n_local(int *n_local, int ncpu, int n){
  
  // put your code here

  
#ifdef DEBUG
  printf("FROM %s: \n",__func__);
#endif

}


void initializeA(float **A, int n, int *n_local, int rank, int ncpu, int argc, char *argv[]){
  
  // optind=1; // hint

   // put your code here
  
#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}


void  create_reference_sol(float **r, int *n_local, int n, int m, int rank, int ncpu, int argc, char *argv[]){

  // put your code here
  
#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}


void parallel_matrix_multiply(float **A, float **r, float **y, int n, int m, int *n_local, int rank, int ncpu){


  // put your code here
  
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
  
  struct { 
    float number; 
    int   index; 
  } localMax, globalMax; 

  // put your code here
  // First find the max index and then use one of the below functions to swap if needed.
  // Note that the row being worked on "might" be the max value,... thus no swapping
  // Needed functions swapLocal and parallelSwap
  
#ifdef DEBUG
  printf("[%d]FROM %s: \n",rank,__func__);
#endif

}

void swapLocal(int global1, int global2, float **A, float **y, int ncpu, int n, int m){


  
  
  // Swap A

  // Swap y


  
#ifdef DEBUG
  printf("[%d]FROM %s: \n",global1%ncpu,__func__);
#endif

}


void parallelSwap(int global1, int global2,float **A, float **y, int ncpu, int n, int m){

  // Swap A and Y

#ifdef DEBUG
  printf("[%d]FROM %s: \n",global1%ncpu,__func__);
#endif

}


void update_rows(float **A,float **y,int n,int m,int *n_local,int rank,int ncpu,int active_row){

  // put your code here
  

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

    // put your code here
  
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
