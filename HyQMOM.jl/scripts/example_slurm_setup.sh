#!/bin/bash
#SBATCH -J hyqmom_setup
#SBATCH -A gts-sbryngelson3
#SBATCH -p cpu-gnr
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=16GB

# ==============================================================================
# HyQMOM Setup for Headless HPC
# ==============================================================================
# This script sets up HyQMOM on a headless cluster (no X11/OpenGL).
# It removes visualization packages and configures MPI.
#
# Usage:
#   1. Edit #SBATCH directives above for your cluster
#   2. sbatch scripts/example_slurm_setup.sh
#   3. Check output: cat slurm-<jobid>.out
# ==============================================================================

echo "Starting HyQMOM setup for headless HPC..."
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo ""

# Load modules (CUSTOMIZE FOR YOUR CLUSTER)
module load julia/1.11
module load openmpi/4.1.5

# Use exec-capable depot (CRITICAL - not HOME if it's noexec!)
export JULIA_DEPOT_PATH="/scratch/$USER/.julia"
mkdir -p $JULIA_DEPOT_PATH

echo "Julia depot: $JULIA_DEPOT_PATH"
echo "Julia version: $(julia --version)"
echo "MPI module: $(module list 2>&1 | grep -i mpi)"
echo ""

cd HyQMOM.jl

# Step 1: Remove visualization packages (headless incompatible)
echo "="^70
echo "STEP 1: Removing visualization packages..."
echo "="^70
julia --project=. scripts/setup_headless.jl

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to remove visualization packages!"
    exit 1
fi

# Step 2: Install remaining packages (including JLD2)
echo ""
echo "="^70
echo "STEP 2: Installing core packages..."
echo "="^70
julia --project=. -e 'using Pkg; Pkg.instantiate()'

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install packages!"
    exit 1
fi

# Step 3: Configure MPI
echo ""
echo "="^70
echo "STEP 3: Configuring MPI..."
echo "="^70
julia --project=. scripts/setup_mpi.jl

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to configure MPI!"
    exit 1
fi

# Step 4: Test
echo ""
echo "="^70
echo "STEP 4: Testing HyQMOM..."
echo "="^70
julia --project=. -e 'using HyQMOM; println("✓ HyQMOM loaded successfully!")'

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to load HyQMOM!"
    exit 1
fi

# Summary
echo ""
echo "="^70
echo "✓ SETUP COMPLETE!"
echo "="^70
echo "Julia depot: $JULIA_DEPOT_PATH"
echo "GLMakie: REMOVED (headless)"
echo "JLD2: INSTALLED (snapshots)"
echo "MPI: CONFIGURED (system)"
echo ""
echo "Next: Submit simulation jobs using scripts/example_slurm_run.sh"
echo "="^70

