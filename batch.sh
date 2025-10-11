#!/bin/bash
#SBATCH -J hyqmom
#SBATCH -A gts-sbryngelson3
#SBATCH -C amd
#SBATCH -N 2                      # number of nodes
#SBATCH --ntasks-per-node=64       # 8 MPI ranks per node (total 16)
#SBATCH --time=00:30:00
#SBATCH --exclusive

module load julia/1.11

# Optional: 1 thread per MPI rank (tune as needed)
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Run: mpirun will respect the Slurm allocation
srun --mpi=pmi2 julia --project=. src/main.jl --no-matplotlib --Np 600 --tmax 0.03