// Program to create Matrix elements for prog05

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

int main(int argc, char *argv[]);

/**********************************************************/
int main(int argc, char *argv[]) {

  int n, row, col, lastColumn, Val;

  n=atoi(argv[1]); // "n" for matrix size nxn
  // argv[2]   create different kinds matrix elements
  // argv[3]   outputfilename

  // example:
  // ./createA 12 f flippedDia.dat

  FILE *fp;
  fp = fopen(argv[3], "w");
  if (fp == NULL){
    fprintf(stderr,"Error while trying to open file %s\n", argv[3]);
    exit(3);
  }

  // First line has n value with a return 
  fprintf(fp,"%d\n",n);  

  
  switch(argv[2][0]) {  // types f, t
  case 'f': // diagonal fliped right-to-left
    /*  1 0 0          0 0 1
	0 2 0  becomes 0 2 0
	0 0 3          3 0 0 */
    lastColumn=n-1;
    for(row=0; row<n; ++row){
      for(col=0; col<n; ++col){
	if(col==lastColumn){
	  // notice no return and a space
	  fprintf(fp,"%.1f ",(float)(row+1));
	  --lastColumn;
	}else{
	  fprintf(fp,"0.0 ");
	}
      }
    }
    fprintf(fp,"\n");  // must end with a return
    fclose(fp);
    break;  // break out of switch stmt
  case 't': // diagonal fliped top-to-bottom
    /*  1 0 0          0 0 3
	0 2 0  becomes 0 2 0
	0 0 3          1 0 0 */
    lastColumn=n-1;  Val=n;
    for(row=0; row<n; ++row){
      for(col=0; col<n; ++col){
	if(col==lastColumn){
	  // notice no return and a space
	  fprintf(fp,"%.1f ",(float)(Val));
	  --lastColumn;
	  --Val;
	}else{
	  fprintf(fp,"0.0 ");
	}
      }
    }
    fprintf(fp,"\n");  // must end with a return
    fclose(fp);    
    break;  // break out of switch stmt
  case 'p': // test for pivot  .. always piviot with bottom row
    /*  1 n+1 1 1 1 1 1..... 1
	2 2 n+1 2 2 2 2 .....2
	3 3 3 n+1 3 3 3 .....3
	:
	n+1 n n n n n n .....n */

    for(row=0; row<(n-1); ++row){
      for(col=0; col<n; ++col){
	if(col==(row+1)){
	  fprintf(fp,"%.1f ",(float)(n+1));
	}else{
	  fprintf(fp,"%.1f ",(float)(row+1));
	}
      }
    }
    row=n-1;
    fprintf(fp,"%.1f ",(float)(n+1));
    for(col=1; col<n; ++col)
      fprintf(fp,"%.1f ",(float)(row+1));

    fprintf(fp,"\n");  // must end with a return
    fclose(fp);
    break;  // break out of switch stmt
  }
  /*
  0.00 0.01 0.02 ...........0.99
  1.00 1.01 1.02 ...........1.99
    :
  99.00 99.01 99.02 ...........99.99
  */
}
