
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



// step 1 obatin n (and m if grad student)
// step 2 calculate n_local for each rank
//        might be helpful for all ranks to know how
//        many rows a give rank has
//        n_local might be an vector....
//        maybe scatter / gather ?

// step 3 use n_local and m to create local memory
// A, r, y, x, difsq, and if rank 0 epsilon 


// step 4 read in the data for A

// step 5 create refference solution

// step 6 Need function to calculate local rhs (y) in parallel
// A*r=y

// Need a funtion to do pivoting in paralel for a given row index
//       This function could use another function
//       That swaps any two rows in parallel
//       Note: the 2 rows could be in the same rank memory


// Need a function to broadcast a row starting at a specific element

// Need a function to update a chosen row and rhs values
//        given a broadcast row

// step 7
// need a function to do gaussian Elimination
//        *** uses above functions ***

// need a function to find next solution element and
//        then update matrix and rhs
//        in parallel

// step 8
// need back substitution method using above functions

// step 9
// function to calcuate error in paralel

// print results

// step 10 *********
// create a testing plan to test varrious inputs knowing
// what the output should be


// use data to create report and reread project to make
// sure nothing was missing






