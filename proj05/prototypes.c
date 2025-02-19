int main(int argc, char *argv[]);
void usage(void);
void get_n_m(int argc, char *argv[], int nAndm[2]);
void define_n_local(int *n_local, int ncpu, int n);
void initializeA(float **A, int n, int *n_local, int rank, int ncpu, int argc, char *argv[]);
void create_referece_sol(float **r, int *n_local, int m, int rank, int ncpu, int argc, char *argv[]);
void parallel_matrix_muliply(float **A, float **r, float **y, int n, int m, int *n_local, int rank, int ncpu);
void pivot(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu, int active_row);
void swap_max_index(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu, int active_row);
void swapLocal(int global1, int global2,float **A, float **y, int ncpu, int n, int m);
void parallelSwap(int global1, int global2,float **A, float **y, int ncpu, int n, int m);
void gaussian_elemination(float **A, float **y, int n, int m, int *n_local, int rank, int ncpu);
void update_rows(float **A,float **y,int n,int m,int *n_local,int rank,int ncpu,int active_row);
void back_substitution(float **A, float **y, float **x, int n, int m, int rank, int ncpu);
void rms_error(float **x, float **r, int *n_local, int n, int m, int rank);
