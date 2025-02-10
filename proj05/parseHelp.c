
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>

#define DEBUG 

int main(int argc, char *argv[]);
void get_n_m(int argc, char *argv[], int nAndm[2]);

/**********************************************************/
int main(int argc, char *argv[]) {

  // *******************************************//
  // step 1 obatin n (and m if grad student)
  //  
  int nAndm[2], n, m;
  get_n_m(argc, argv, nAndm);
  n=nAndm[0];         
  m=nAndm[1];         
  printf("n=%d m=%d\n",n,m);

}

void get_n_m(int argc, char *argv[], int nAndm[2]){
  int c, n=-1, m=10;
  FILE *fp;


  //  optind=1;
  int stuff_goes_here;
  while(1){
    
    c = getopt(argc, argv, "r:d:f:R:F:");  
    
    if( c == -1) {  
      break; // nothing left to parse in the command line
    }
    
    switch(c) { 
    case 'r': n=atoi(optarg); 
      break;  // break out of switch stmt
    case 'd': n=atoi(optarg);
      break;  // break out of switch stmt
    case 'f': 
      stuff_goes_here=0;
      break;  // break out of switch stmt
    case 'R':
      stuff_goes_here=0;
      break;  // break out of switch stmt
    case 'F': 
      stuff_goes_here=0;
      break;  // break out of switch stmt
    default: fprintf(stderr,"This flag DNE: %c\n", optopt);
      exit(1);
    }
  }
  
  if (n<1){
    fprintf(stderr," n was not defined .... code stoped\n");
    exit(2);
  }
  
  nAndm[0]=n;
  nAndm[1]=m;

#ifdef DEBUG
  printf("FROM %s: n=%d m=%d\n",__func__,n,m);
#endif

}
