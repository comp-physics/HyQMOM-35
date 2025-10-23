#!/bin/bash
#SBATCH -J hyqmom
#SBATCH -A gts-sbryngelson3
#SBATCH -p cpu-small
#SBATCH -N 2                      # number of nodes
#SBATCH --ntasks-per-node=8       # 8 MPI ranks per node (total 16)
#SBATCH --time=00:30:00

# --exclusive

module load julia/1.11

# Optional: 1 thread per MPI rank (tune as needed)
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1

export JULIA_DEPOT_PATH="$HOME/.julia"

# Change to Julia package directory
cd HyQMOM.jl

srun --nodes=1 --ntasks=1 julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# 2) Disable per-rank auto precompile spam during the big run
export JULIA_PKG_PRECOMPILE_AUTO=0

# Run: mpirun will respect the Slurm allocation
srun --mpi=pmi2 julia --project=. examples/run_3d_custom_jets.jl --no-viz true --Nx 50 --Ny 50 --Nz 50 --tmax 0.01 --Ma 1.2 --Kn 0.1  --snapshot-interval 1
