//
// cs4900
// Project 04
// Shapiy Sagiev
// w1140sxs
// Due 22 Feb 2024, accepted up to Feb 23 10pm
// System = fry
// Compiler syntax = ./compile.sh proj04b.c
// Job Control File = proj04b.sbatch
// Additional File  = NA
// Results file     = proj04b.txt
//

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

float dot(float *v1, float *v2, int N) {
  int i;
  float sum = 0.0;
  for (i = 0; i < N; i++) {
    sum += v1[i] * v2[i];
  }
  return sum;
}

void multiply(float *v1, float *v2, int Nr, int L1, float **matrix) {
  for (int i = 0; i < Nr; i++) {
    for (int j = 0; j < L1; j++) {
      matrix[i][j] = v1[i] * v2[j];
    }
  }
}

int main(int argc, char *argv[]) {
  // Inputs
  int option =
      atoi(argv[1]); // L1 determines whether to run vector math or matrix math
  int A = 1;         // Variable hardset because does not change each run
  int B, Ap, Bp;     // Variables not set because they change with each loop

  // Standard MPI initialization calls
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int rank, ncpu, name_len;
  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
  MPI_Get_processor_name(processor_name, &name_len);

  if (option > 0) // If First arguement is negative, do only vector math...
  {
    if (rank == 0) { printf("\n===========================VECTOR MATH===========================\n"); }
    int flag = 1;
    for (int p = 1; p <= 9; p++) // 9 different values of p
    {
      int L1 = pow(10, p); // Assign L1 to 10^p
      if (option == 1) {
        if (rank == 0 && flag) {
          printf("==============OPTION 1 RUN==============");
          flag = 0;
        }
        B = -1;
        Ap = 45;
        Bp = 45;
      } else if (option == 2) {
        if (rank == 0 && flag) {
          printf("==============OPTION 2 RUN==============");
          flag = 0;
        }
        B = (int)(L1 / 4);
        Ap = 22;
        Bp = 90;
      } else if (option == 3) {
        if (rank == 0 && flag) {
          printf("==============OPTION 3 RUN==============");
          flag = 0;
        }
        B = (int)(2.2 * L1);
        Ap = 75;
        Bp = 135;
      } else {
        if (rank == 0 && flag) {
          printf("==============CUSTOM OPTION RUN==============");
          flag = 0;
        }
        B = atoi(argv[2]);
        Ap = atoi(argv[3]);
        Bp = atoi(argv[4]);
      }
      // Constants
      float pi, dx, norm, a, b, Xi;
      pi = acos(-1.0);
      dx = 2.0 * pi / L1;
      norm = sqrt(2.0 / L1);
      // Phase offset
      a = pi * Ap / 180;
      b = pi * Bp / 180;

      // Calculates which processor has what part of the global elements
      int Nr, R;
      float xlast;
      Nr = (int)(L1 / ncpu); // All processors have at least this many rows
      R = (L1 % ncpu);       // Remainder of above calculation
      if (rank < R) {
        ++Nr; // First R processors have one more element
        xlast = -dx / 2.0 + rank * Nr * dx;
      } else {
        xlast = -dx / 2.0 + R * (Nr + 1) * dx + (rank - R) * Nr * dx;
      }

      // Initialize the elements of v1 and v2 on this partition
      float *v1 = (float *)malloc(sizeof(float) * Nr);
      assert(v1 != NULL);
      float *v2 = (float *)malloc(sizeof(float) * Nr);
      assert(v2 != NULL);
      // Vectors 1, 2 Filled
      Xi = xlast;
      for (int i = 0; i < Nr; i++) {
        Xi = Xi + dx;
        v1[i] = norm * cos(A * Xi + a);
        v2[i] = norm * sin(B * Xi + b);
      }

      // Dot Product
      float local_sum;
      local_sum = dot(v1, v2, Nr);
      // Reduce all of the local sums into the global sum
      float global_sum;
      MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

      if (rank == 0) {
        printf("\nStarting to process vector dot sum with L1=%i, A=%i, B=%i, "
               "Ap=%i, Bp=%i",
               L1, A, B, Ap, Bp);
        printf("\nL1=10^%i: ", p);
        printf("\nThe vector lengths are %d elements long", Nr);
        printf("\nThe dot product is %f\n", global_sum);
      }
      // Free Memory
      free(v1); // Both Vector and Matrix math use v1
      free(v2);
    }
  } else // Matrix Math
  {
    if (rank == 0) { printf("\n===========================MATRIX MATH===========================\n"); }
    int flag = 1;
    for (int p = 1; p <= 100; p += 99) // 9 different values of p
    {
      // Changing Variables
      int L1 = 1000 * p;
      if (option == -1) {
        if (rank == 0 && flag) {
          printf("==============OPTION 1 RUN==============");
          flag = 0;
        }
        B = 1;
        Ap = 135;
        Bp = 135;
      } else if (option == -2) {
        if (rank == 0 && flag) {
          printf("==============OPTION 2 RUN==============");
          flag = 0;
        }
        B = (int)(2.4 * L1);
        Ap = 45;
        Bp = 135;
      } else {
        if (rank == 0 && flag) {
          printf("==============CUSTOM OPTION RUN==============");
          flag = 0;
        }
        B = atoi(argv[2]);
        Ap = atoi(argv[3]);
        Bp = atoi(argv[4]);
      }
      // Constants
      float pi, dx, norm, a, b;
      pi = acos(-1.0);
      norm = sqrt(2.0 / L1);
      dx = 2.0 * pi / L1; // dy=dx because matrix has same L1=L2, no difference between vector sizes, so no need for dy
      // Phase offset
      a = pi * Ap / 180;
      b = pi * Bp / 180;
      // Calculates which processor has what part of the global elements
      // Block row Partitioned
      int Nr, R, row;
      float xlast, Xi;
      Nr = (int)(L1 / ncpu); // All processors have at least this many rows
      R = (L1 % ncpu);       // Remainder of above calculation
      if (rank < R) {
        ++Nr;                // First R processors have one more element
        xlast = -dx / 2.0 + rank * Nr * dx;
      } else {
        xlast = -dx / 2.0 + R * (Nr + 1) * dx + (rank - R) * Nr * dx;
      }

      // Initialize the matrix rows on this partition
      // and the global vector
      int i, j;
      float **matrix = (float **)malloc(Nr * sizeof(float *));
      assert(matrix != NULL);
      for (i = 0; i < Nr; i++) {
        matrix[i] = (float *)malloc(sizeof(float) * L1);
        assert(matrix[i] != NULL);
      }

      // Initialize the elements of v1 and v2 on this partition
      float *v1 = (float *)malloc(sizeof(float) * Nr);
      assert(v1 != NULL);
      float *v3 = (float *)malloc(sizeof(float) * L1);
      assert(v3 != NULL);

      // Vectors 1, 3 Filled
      Xi = xlast;
      for (i = 0; i < L1; i++) {
        if (i < Nr) // Since v1 is length of Nr, fill vector as long as i < Nr
        {
            v1[i] = norm * cos(A * Xi + a);
        }
        v3[i] = norm * cos(B * Xi + b);
      }
      // Multiplying vectors to fill matrix
      multiply(v1, v3, Nr, L1, matrix);
      // Each rank will go through its diags
      float local_diag_sum = 0;
      for (j = 0; j < Nr; j++) {
        local_diag_sum += matrix[j][j + (Nr * rank)]; // Nr * rank adjusts for starting row
      }
      float global_diag_sum;
      MPI_Reduce(&local_diag_sum, &global_diag_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      // Sum up all elements
      float local_sum = 0;
      for (i = 0; i < Nr; i++) {
        for (j = 0; j < L1; j++) {
          local_sum += matrix[i][j];
        }
      }
      local_sum -= local_diag_sum; // Subtract diagonal row values
      float global_sum;
      MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
      if (rank == 0) {
        printf("\nMatrix Results");
        printf("\nMatrix Size: [%i x %i]", L1, L1);
        printf("\nPartitions handle at least %i Rows.", Nr);
        printf("\nThe diagonal sum is: %f\n", global_diag_sum);
        printf("Matrix total sum (w/o diagonal sum) is: %f", global_sum);
      }
      // Free memory
      free(v1); // Both Vector and Matrix math use v1
      free(v3);
      for (i = 0; i < Nr; i++) {
        free(matrix[i]); // Free each row of the matrix
      }
      free(matrix); // Free the matrix array itself
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
