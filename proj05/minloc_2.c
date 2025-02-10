#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#define  N   100

struct { 
  float number; 
  int   index; 
} localMin, globalMin; 



int main(int argc, char** argv) {
  
  float vector[N];        /* local array of random numbers */ 
  int minRank, minIndex; 
  float minNumber; 
    
  int rank, ncpu, root=0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

  /*
  // time(NULL) returns a different number each time called
  // -10*rank makes sure each rank has a different seed
  printf("[%d] %d\n",rank,time(NULL));
  printf("[%d] %d\n",rank,time(NULL)-10*rank);
  */

  srand(time(NULL)+10*rank);  // each rank gets a different seed
  
    // Put some data into the vector in
  // order to test the rest of the code
  int i;
  for (i=0; i<N; ++i){
    vector[i]=(float)( 1000000.0*( ((rand()%2147483647)/2147483647.0)-0.5) );
  }
  //vector[i] = (float)( (rand()%(20*N)-10*N) );
    
  

  // *************************************************//
  //                                                  //
  // Start main idea here                             //
  //                                                  //


  /* local minloc */ 
  localMin.number = vector[0]; 
  localMin.index = 0; 
  for (i=1; i < N; i++) 
    if (localMin.number > vector[i]) { 
      localMin.number = vector[i]; 
      localMin.index = i; 
    } 

  printf("[%d] minval=%f localIndex=%d\n",rank,localMin.number,localMin.index);
  
  /* local to global index transformation */ 
  localMin.index = rank + localMin.index*ncpu; 
  MPI_Reduce(&localMin, &globalMin, 1, MPI_FLOAT_INT, MPI_MINLOC, root, MPI_COMM_WORLD ); 
  /* At this point, the answer resides on process root 
   */ 

  //                                                  //
  // End main idea here                               //
  //                                                  //
  // *************************************************//

    
  if (rank == root) { 
    minNumber = globalMin.number; 
    minRank = globalMin.index%ncpu; 
    minIndex = (globalMin.index-minRank)/ncpu; 
    printf("[%d] Result: minval=%f globalIndex= %d minrank=%d localIndex=%d\n",rank,minNumber,globalMin.index,minRank,minIndex);
  }

  MPI_Finalize();
 
}
