#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj04
#SBATCH --output=proj04.txt
#SBATCH --time=5:00

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0

## compile your code
mpicc proj04a.c -o proj04.mpi -lm

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Make code executable
chmod u+x proj04.mpi

# Run on a a specified number of cores
# mpirun -np 4 ./proj04.mpi 1000
./proj04.mpi 8 10 2 3 0 0
