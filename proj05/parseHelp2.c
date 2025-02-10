
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


//#define DEBUG 

int main(int argc, char *argv[]);
void get_n_m(int argc, char *argv[], int nAndm[2]);

/**********************************************************/
int main(int argc, char *argv[]) {

  // *******************************************//
  // step 1 obatin n (and m if grad student)
  //  
  int nAndm[2], n, m;
  printf("Output for the first call\n");
  get_n_m(argc, argv, nAndm);

  printf("\n\nOutput for 2nd call.. without setting optind=1\n");
  // Call it again without setting optind=1;
  get_n_m(argc, argv, nAndm);

  printf("\n\nOutput for 3rd call.. but first set optind=1\n");
  // Call it again with optind=1;
  optind=1;
  get_n_m(argc, argv, nAndm);

  
}

void get_n_m(int argc, char *argv[], int nAndm[2]){
  int c, n=-1, m=10;
  FILE *fp;

  nAndm[0]=n;         
  nAndm[1]=m;         


  //  optind=1;
  int stuff_goes_here;
  while(1){
    
    c = getopt(argc, argv, "a:b:c:Q:W");  
    
    if( c == -1) {  
      break; // nothing left to parse in the command line
    }
    
    switch(c) { 
    case 'a':
      printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg); 
      break;  // break out of switch stmt
    case 'b':
      printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg); 
      break;  // break out of switch stmt
    case 'c':
      printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg); 
      break;  // break out of switch stmt
    case 'Q':
      printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg); 
      break;  // break out of switch stmt
    case 'W':
      printf("Entry number %d was -%c and it had a value of %s\n",optind,c,optarg); 
      break;  // break out of switch stmt
    default: fprintf(stderr,"This flag DNE: %c\n", optopt);
      //exit(1);
    }
  }
  printf("get_n_m done ...... \n");
  
#ifdef DEBUG
  printf("FROM %s: n=%d m=%d\n",__func__,n,m);
#endif

}
