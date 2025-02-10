#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=proj04
#SBATCH --output=proj04b.txt
#SBATCH --time=30:00

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0

## compile your code
mpicc proj04b.c -o proj04b.mpi -lm

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Make code executable
chmod u+x proj04b.mpi

# Run on a a specified number of cores
for arg in {1..3}; do
	mpirun -np 16 ./proj04b.mpi $arg >> proj04b.txt
done
for arg in $(seq -1 -1 -2); do
	mpirun -np 16 ./proj04b.mpi $arg >> proj04b.txt
done
