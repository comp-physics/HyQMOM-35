#!/bin/bash
#SBATCH -J hyqmom
#SBATCH -A gts-sbryngelson3
#SBATCH -p cpu-small,cpu-medium,cpu-large
#SBATCH -N 6                       # number of nodes
#SBATCH --ntasks-per-node=16       # 8 MPI ranks per node (total 16)
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=5GB

#--constraint=graniterapids

module load julia/1.11

export HYQMOM_SKIP_PLOTTING=true
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_PKG_PRECOMPILE=0

# Optional: 1 thread per MPI rank (tune as needed)
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1

export JULIA_DEPOT_PATH="$HOME/.julia"

# Change to Julia package directory
cd HyQMOM.jl

# srun --nodes=1 --ntasks=1 julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
srun --nodes=1 --ntasks=1 julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 2) Disable per-rank auto precompile spam during the big run
export JULIA_PKG_PRECOMPILE_AUTO=0

# Run: mpirun will respect the Slurm allocation
srun --mpi=pmi2 julia --project=. examples/run_3d_custom_jets.jl --Nx 200 --Ny 200 --Nz 200 --tmax 0.02 --Ma 1.0 --Kn 0.1  --snapshot-interval 1 --no-viz true
