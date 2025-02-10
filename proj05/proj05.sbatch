#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=3
#SBATCH --job-name=proj05
#SBATCH --output=proj05.txt
#SBATCH --time=10:00

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0

## compile your code
mpicc proj05.c -g -o proj05.mpi -lm

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Make code executable
chmod u+x proj05.mpi

# Run on a a specified number of cores
mpirun -np 3 ./proj05.mpi -r 10
#mpiexec -np 3 valgrind --tool=memcheck --leak-check=full ./proj05.mpi -r 10

