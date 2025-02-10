//
// cs4900
// Project 03
// Dr. John Nehrbass
// w006jwn
// Due 08 Feb 2024, accepted up to Feb 09 10pm
// System = fry
// Compiler syntax = ./compile.sh proj03
// Job Control File = proj03.sbatch
// Additional File  = NA
// Results file     = proj03.txt
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

//#define DEBUG

// DO NOT CHANGE THIS FUNCTION!!!!!
// function to return a random variable between [0,1)
double randVar(){

  // Notes:
  // 2147483647 is the largest possible value of rand()
  // printf("max=%llu\n",RAND_MAX);

  return (double)((rand()%2000000000)/2000000000.0);
}


/* function to return M and P*/
unsigned long* calMandP(int N ) {

  //Save the results in an integer array
  static unsigned long MnP[2];

  // Constant for small square 
  // Reciprocal of sqrt(2) or 1/sqrt(2)
  double rsq2=0.70710678118654752440;

  // initialize counters each time the function is called
  unsigned long M=0, P=0;
  
  // random variables
  double x, y, rs;
  int randomPnts;

   // Pick N random (x,y) points
  for(randomPnts = 0; randomPnts < N; ++randomPnts) {

    x=randVar();
    y=randVar();

    // radius squared
    rs=x*x+y*y;

#ifdef DEBUG
    // test random values in function
    printf("%d x=%f y=%f rs=%f\n",randomPnts,x,y,rs);
#endif 

    if (rs<1){                   // or should it be rs<=1
      // Inside circle
      ++M;
      if (x<rsq2 && y<rsq2){     // or should it be <=
				// inside inner square
				++P;
      }
    }
  }

  MnP[0]=M;  
  MnP[1]=P;  
  return MnP;
}

// Global Variables
int count = 0; // For how many times looped
unsigned long total[2]; // Total points "thrown"
unsigned long mTotal = 0; // Totals of M and P
unsigned long pTotal = 0;
int main (int argc, char *argv[]){
  clock_t begin = clock();
	int rank, ncpu;
  char my_cpu_name[1024];
  int my_name_length;
	
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
  MPI_Get_processor_name(my_cpu_name, &my_name_length);
 
	float epsilon = 0.001; // Epsilon Value
	bool repeat = true;
	while (repeat) {
		// Rank 0 CPU is the "Main" CPU which will track repeats and perform calculations
		if (rank == 0) {
			count++;
			printf("\nLoop #%i", count);
		}
  	// For convenience I passed "N" on the command line
  	// Not checking for valid input
  	char *a = argv[1];
  	int N = atoi(a);
		// Use current time as seed for random generator
  	srand(rank+count);
  
  	// The Max value of unsigned long is 18,446,744,073,709,551,615
  	unsigned long* MnP;   // Store in vector M=MnP[0] and P=MnP[1]
  	// Return M and P for this case
  	MnP=calMandP(N);
		MPI_Barrier(MPI_COMM_WORLD);
		// Once all points are calculated...
		if (rank != 0) {
			MPI_Send(MnP,2,MPI_UNSIGNED_LONG,0,0,MPI_COMM_WORLD);
		} else {
			total[0] += MnP[0]; // Add rank 0 MnP to total
			total[1] += MnP[1];
			for (int i=1; i<4; i++) { // Loop thru CPUs 1-3 and recieve their MnP values
				MPI_Recv(MnP,2,MPI_UNSIGNED_LONG,i,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				total[0] += MnP[0];
				total[1] += MnP[1];
			}
			// Results (Calculated by rank 0 CPU)
      printf("\nM=%llu P=%llu\n",total[0],total[1]);
      printf("use(pi=4M/N) pi~%f\n",4.0*total[0]/(4*N*count));
      printf("use(pi=2M/P) pi~%f\n",2.0*total[0]/total[1]);

			if ((((4.0*mTotal/(4*N*(count-1)))-(4.0*total[0]/(4*N*count)) >= epsilon) &&
					 (((2.0*total[0]/total[1])-(2.0*mTotal/pTotal)) <= epsilon)) && (count != 1)) {
				printf("\nResults\nN: %i\n", N);
				printf("Iterations: %i\n", count);
				printf("M Total: %llu\n", mTotal);
				printf("N Total: %llu\n", pTotal);
				printf("Epsilon: %f\n", epsilon);
				
				repeat=false;
			}
			// Former results saved for comparison
      mTotal = total[0];
      pTotal = total[1];
		}
	}
	
	clock_t end = clock();
	double time = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time Spent: %f\n", time);
	MPI_Finalize();
}
