#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=minloc_2
#SBATCH --output=minloc_2.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0

## compile your code
mpicc minloc_2.c -o minloc_2.mpi -lm

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Make code executable
chmod u+x minloc_2.mpi

# Run on a single Node
echo "run one"
mpirun -np 3 ./minloc_2.mpi
sleep 1
echo ; echo; echo;
echo "run two"
mpirun -np 3 ./minloc_2.mpi
sleep 1
echo ; echo; echo;
echo "run three"
mpirun -np 3 ./minloc_2.mpi
sleep 1
echo ; echo; echo;
echo "run four"
mpirun -np 3 ./minloc_2.mpi
sleep 1
echo ; echo; echo;
echo "run five"
mpirun -np 16 ./minloc_2.mpi
