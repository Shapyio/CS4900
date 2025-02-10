#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=minloc
#SBATCH --output=minloc.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0

## compile your code
mpicc minloc.c -o minloc.mpi -lm

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Make code executable
chmod u+x minloc.mpi

# Run on a single Node 
mpirun -np 3 ./minloc.mpi
sleep 1
mpirun -np 3 ./minloc.mpi
sleep 1
mpirun -np 3 ./minloc.mpi
sleep 1
mpirun -np 3 ./minloc.mpi
sleep 1
mpirun -np 16 ./minloc.mpi
