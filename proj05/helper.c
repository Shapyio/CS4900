// helper functions
void printA(float **A, int n, int rank, int cpu);
void printr(float **r, int n, int m, int rank, int cpu);
void printx(float **x, int n, int m, int rank, int cpu);
void printy(float **y, int n, int m, int rank, int cpu);

// functions

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
