#!/bin/bash
#SBATCH -J hyqmom_Ma70_Kn10
#SBATCH -A gts-sbryngelson3
#SBATCH -p cpu-gnr
#SBATCH -N 8
#SBATCH --ntasks-per-node=96
#SBATCH --time=8:00:00
#SBATCH --mem=128GB
#SBATCH --constraint=graniterapids
#SBATCH --qos=embers

# ==============================================================================
# HyQMOM Production Run on HPC
# ==============================================================================
# This script runs a production HyQMOM simulation with MPI parallelization.
#
# Prerequisites:
#   - Setup completed: scripts/example_slurm_setup.sh
#   - LocalPreferences.toml exists with MPI configuration
#
# Usage:
#   1. Edit #SBATCH directives and parameters below
#   2. sbatch scripts/example_slurm_run.sh
#   3. Monitor: tail -f slurm-<jobid>.out
# ==============================================================================

echo "Starting HyQMOM simulation..."
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks: $SLURM_NTASKS"
echo ""

# Load modules (SAME as setup script!)
module load julia/1.11
module load openmpi/4.1.5

# Use SAME depot as setup
export JULIA_DEPOT_PATH="/scratch/$USER/.julia"

# Disable plotting and precompilation
export HYQMOM_SKIP_PLOTTING=true
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1

echo "Julia depot: $JULIA_DEPOT_PATH"
echo "HYQMOM_SKIP_PLOTTING: $HYQMOM_SKIP_PLOTTING"
echo ""

cd HyQMOM.jl

# Check MPI configuration
if [ ! -f "LocalPreferences.toml" ]; then
    echo "WARNING: LocalPreferences.toml not found!"
    echo "Run setup script first: sbatch scripts/example_slurm_setup.sh"
    exit 1
fi

echo "MPI Configuration:"
cat LocalPreferences.toml
echo ""

# ============================================================================
# SIMULATION PARAMETERS - EDIT THESE
# ============================================================================
NX=80           # Grid points in x
NY=80           # Grid points in y
NZ=80           # Grid points in z
TMAX=0.0025     # Simulation end time
MA=70.0         # Mach number
KN=10.0         # Knudsen number
CFL=0.3         # CFL number (default: 0.3)
SNAPSHOT=50     # Snapshot interval (steps)
CONFIG="crossing"  # Jet configuration: crossing, crossing2D, triple-jet, quad-jet

echo "Simulation Parameters:"
echo "  Grid: ${NX}x${NY}x${NZ}"
echo "  tmax: $TMAX"
echo "  Ma: $MA"
echo "  Kn: $KN"
echo "  CFL: $CFL"
echo "  Snapshot interval: $SNAPSHOT steps"
echo "  Configuration: $CONFIG"
echo ""

# Run simulation with MPI
echo "Starting simulation with MPI..."
echo "="^70

srun --mpi=pmix julia --project=. examples/run_3d_custom_jets.jl \
  --Nx $NX --Ny $NY --Nz $NZ \
  --tmax $TMAX \
  --Ma $MA \
  --Kn $KN \
  --CFL $CFL \
  --snapshot-interval $SNAPSHOT \
  --config $CONFIG \
  --no-viz true

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

