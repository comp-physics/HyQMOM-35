#!/bin/bash
#SBATCH -J hyqmom
#SBATCH -A gts-sbryngelson3
#SBATCH -p cpu-gnr
#SBATCH -N 8                       # number of nodes
#SBATCH --ntasks-per-node=96       # 8 MPI ranks per node (total 16)
#SBATCH --time=8:00:00
#SBATCH --mem=128GB
#SBATCH --constraint=graniterapids
#SBATCH --qos=embers

module load julia/1.11
module load openmpi/4.1.5

# Use SAME depot as setup
export JULIA_DEPOT_PATH="$HOME/scratch/.julia"

# Disable plotting and precompilation
export HYQMOM_SKIP_PLOTTING=true
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMPI_MCA_btl=^openib
export UCX_WARN_UNUSED_ENV_VARS=n


echo "Julia depot: $JULIA_DEPOT_PATH"
echo "HYQMOM_SKIP_PLOTTING: $HYQMOM_SKIP_PLOTTING"
echo ""

# Ensure we're in HyQMOM.jl directory (script should be run from there)
if [ ! -f "Project.toml" ]; then
    echo "ERROR: Project.toml not found! Must run from HyQMOM.jl directory:"
    echo "  cd HyQMOM.jl"
    echo "  sbatch slurm/hyqmom_base.sbatch"
    exit 1
fi

# Check MPI configuration
if [ ! -f "LocalPreferences.toml" ]; then
    echo "WARNING: LocalPreferences.toml not found!"
    echo "Run setup script first: sbatch slurm/setup_hyqmom.sbatch"
    exit 1
fi

echo "MPI Configuration:"
cat LocalPreferences.toml
echo ""



# Run simulation with MPI
srun --mpi=pmix julia --project=. examples/run_3d_custom_jets.jl --Nx 80 --Ny 80 --Nz 80 --tmax 0.0025 --Ma 70.0 --Kn 10.0 --snapshot-interval 50 --no-viz
EXIT_CODE=$?

echo ""
echo "="^70
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ SIMULATION COMPLETE!"
    echo ""
    echo "Output files:"
    ls -lh snapshots_*.jld2 2>/dev/null || echo "  (no snapshot files found)"
    echo ""
    echo "To visualize results:"
    echo "  1. Copy to local: scp $USER@cluster:$(pwd)/snapshots_*.jld2 ."
    echo "  2. Visualize: julia --project=. visualize_jld2.jl snapshots_*.jld2"
else
    echo "✗ SIMULATION FAILED (exit code: $EXIT_CODE)"
    echo "Check output above for errors"
fi
echo "="^70
echo "Job ended: $(date)"