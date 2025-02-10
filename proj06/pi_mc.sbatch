#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=proj06
#SBATCH --output=proj06.txt
#SBATCH --time=10:00

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0

## Compile your code
mpicc pi_mc.c Random.c -o pi_mc.mpi -lm -fopenmp

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Make code executable
chmod u+x pi_mc.mpi

# Run on a a specified number of cores
mpirun -np 1 ./pi_mc.mpi 1000000 1
