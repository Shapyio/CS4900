//
// cs4900
// Project 04
// Shapiy Sagiev
// w1140sxs
// Due 22 Feb 2024, accepted up to Feb 23 10pm
// System = fry
// Compiler syntax = ./compile.sh proj04a
// Job Control File = proj04.sbatch
// Additional File  = NA
// Results file     = proj04.txt
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

float dot(float *v1, float *v2, int N) {
  int i;
  float sum=0.0;
  for (i = 0; i < N; i++) {
    sum += v1[i]*v2[i];
  }
  return sum;
}

void multiply(float *v1, float *v2, int L1, int L2, float **matrix){
	for (int i=0; i<L2; i++){
		for (int j=0; j<L1; j++){
			matrix[i][j]=v1[j]*v2[i];
		}
	}
}

int main (int argc, char *argv[]){
	// Inputs
  int L1 = atoi(argv[1]);
  int L2 = atoi(argv[2]);
  int A  = atoi(argv[3]);
  int B  = atoi(argv[4]);
  int Ap = atoi(argv[5]);
	int Bp = atoi(argv[6]);

	// Constants
  float dx, pi, norm, Xi, dy, Yi, a, b;
  pi=acos(-1.0);
  dx=2.0*pi/L1;
  norm=sqrt(2.0/L1);
	dy=2.0*pi/L2;
	// Phase offset
  a=pi*Ap/180;
  b=pi*Bp/180;

	// Initialize the elements on this partition
  float *v1 = (float *)malloc(sizeof(float) * L1);
  assert(v1 != NULL);
  float *v2 = (float *)malloc(sizeof(float) * L1);
  assert(v2 != NULL);
	float *v3 = (float *)malloc(sizeof(float) * L2);
	assert(v3 != NULL);
	
	// Vectors 1, 2
  for (int i = 0; i < L1; i++) {
		Xi=dx/2.0+i*dx;
    v1[i]=norm*cos(A*Xi+a);
    v2[i]=norm*sin(B*Xi+b);
  }
	// Vector 3
	for (int j = 0; j < L2; j++){
		Yi=dy/2.0+j*dy;
		norm=sqrt(2.0/L2); // Normalize for L2
		v3[j]=norm*cos(B*Yi+b);
	}

	// Dot Product
  float local_sum;
  local_sum=dot(v1,v2,L1);
	// Multiplication
	// Allocate memory for the matrix
	float **matrix = (float **)malloc(sizeof(float *) * L2);
	assert(matrix != NULL);
	for (int i = 0; i < L2; i++) {
    matrix[i] = (float *)malloc(sizeof(float) * L1);
    assert(matrix[i] != NULL);
	}
	multiply(v1, v3, L1, L2, matrix);
	
	// Print Results
	printf("--------Proj04a.c Results--------\n");
	printf("Dot Product:\n");
	printf("v1 * v2 = (");
	for (int x=0; x<L1; x++){
		printf("%f ", v1[x]);
	}
	printf(") * (");
	for (int y=1; y<L1; y++){
		printf("%f ", v2[y]);
	}
	printf(")\nResult is: %f", local_sum);
	printf("\nV3 = (");
	for (int z=0; z<L2; z++){
		printf("%f ", v3[z]);
	}
	printf(")\n\nMultiplication: \n");
	for (int i=0; i<L2; i++){
		for (int j=0; j<L1; j++){
			printf("%f ", matrix[i][j]);
		}
		printf("\n\n");
	}
	// Clean up Memory
  free(v1);
  free(v2);
	free(v3);
	free(matrix);
}
