#!/usr/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=2
#SBATCH --job-name=parseHelp
#SBATCH --output=parseHelp.txt

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0


## compile your code
mpicc parseHelp.c -o parseHelp.mpi -lm

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Make code executable
chmod u+x parseHelp.mpi

# Run on a a specified number of cores
mpirun -np 2 ./parseHelp.mpi


