#!/bin/bash

#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH --job-name=report
#SBATCH --output=report.txt
#SBATCH --time=15:00

## First make sure all other modules are not loaded
module purge

## Load the appropriate module
module load openmpi/4.0.3-gcc-8.5.0

# Fix your environment
export OMPI_MCA_pml="^ucx"
export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

# Run on a a specified number of cores
mpirun -np 16 ~/../w006jwn/proj05/report.mpi -r 10000
