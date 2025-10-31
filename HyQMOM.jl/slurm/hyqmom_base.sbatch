#!/bin/bash
# ==============================================================================
# HyQMOM Base Template for Parameter Sweeps
# ==============================================================================
# This is a TEMPLATE file used by slurm/sweep.sh to generate actual job scripts.
# 
# To submit a single job, edit the parameters below and run:
#   cd HyQMOM.jl
#   sbatch slurm/hyqmom_base.sbatch
#
# To submit parameter sweeps, use slurm/sweep.sh instead (recommended).
# ==============================================================================

#SBATCH -J hyqmom_TEMPLATE
#SBATCH -A gts-sbryngelson3
#SBATCH -p cpu-gnr
#SBATCH -N 8                       # number of nodes
#SBATCH --ntasks-per-node=96       # 8 MPI ranks per node (total 16)
#SBATCH --time=8:00:00
#SBATCH --mem=128GB
#SBATCH --constraint=graniterapids
#SBATCH --qos=embers

# Ensure we're in HyQMOM.jl directory (check FIRST before doing anything)
if [ ! -f "Project.toml" ]; then
    echo "ERROR: Project.toml not found! Must run from HyQMOM.jl directory:" >&2
    echo "Usage:" >&2
    echo "  cd HyQMOM.jl" >&2
    echo "  sbatch slurm/hyqmom_base.sbatch" >&2
    exit 1
fi

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

# Verify no visualization packages (headless HPC check)
julia --project=. -e '
using Pkg
deps = keys(Pkg.project().dependencies)
bad_packages = ["GLMakie", "FileIO", "ColorSchemes", "LaTeXStrings", "MAT", "HDF5"]
found = filter(p -> p in deps, bad_packages)
if !isempty(found)
    println(stderr, "")
    println(stderr, "="^70)
    println(stderr, "ERROR: Visualization packages found in Project.toml!")
    println(stderr, "="^70)
    println(stderr, "Found: ", join(found, ", "))
    println(stderr, "")
    println(stderr, "These packages will fail to precompile on headless HPC.")
    println(stderr, "Run the setup script to remove them:")
    println(stderr, "  sbatch slurm/setup_hyqmom.sbatch")
    println(stderr, "="^70)
    exit(1)
end
' || exit 1

echo "✓ Headless environment verified"
echo ""

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
printf '=%.0s' {1..70}; echo ""
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
printf '=%.0s' {1..70}; echo ""
echo "Job ended: $(date)"