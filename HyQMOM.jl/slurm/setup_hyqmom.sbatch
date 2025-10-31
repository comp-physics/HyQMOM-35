#!/bin/bash
#SBATCH -J hyqmom_setup
#SBATCH -A gts-sbryngelson3
#SBATCH -p cpu-gnr
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=16GB

# ==============================================================================
# HyQMOM Setup for Headless HPC
# ==============================================================================
# This script sets up HyQMOM on a headless cluster (no X11/OpenGL).
# It removes visualization packages and configures MPI.
#
# Usage:
#   1. Edit #SBATCH directives above for your cluster
#   2. cd HyQMOM.jl
#   3. sbatch slurm/setup_hyqmom.sbatch
#   4. Check output: cat slurm-<jobid>.out
# ==============================================================================

# Ensure we're in HyQMOM.jl directory (check FIRST before doing anything)
if [ ! -f "Project.toml" ]; then
    echo "ERROR: Project.toml not found! Must run from HyQMOM.jl directory:" >&2
    echo "Usage:" >&2
    echo "  cd HyQMOM.jl" >&2
    echo "  sbatch slurm/setup_hyqmom.sbatch" >&2
    exit 1
fi

echo "Starting HyQMOM setup for headless HPC..."
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo ""

# Load modules (CUSTOMIZE FOR YOUR CLUSTER)
module load julia/1.11
module load openmpi/4.1.5

# Use exec-capable depot (CRITICAL - not HOME if it's noexec!)
export JULIA_DEPOT_PATH="$HOME/scratch/.julia"
mkdir -p $JULIA_DEPOT_PATH

# Disable auto-precompile to avoid stale cache errors
export JULIA_PKG_PRECOMPILE_AUTO=0

echo "Julia depot: $JULIA_DEPOT_PATH"
echo "Julia version: $(julia --version)"
echo "MPI module: $(module list 2>&1 | grep -i mpi)"
echo ""

# Step 0: Clear precompiled cache BEFORE doing anything
echo ""
printf '=%.0s' {1..70}; echo ""
echo "STEP 0: Clearing precompiled cache..."
printf '=%.0s' {1..70}; echo ""
JULIA_VERSION=$(julia --version | grep -oP '\d+\.\d+' | head -1)
CACHE_DIR="$JULIA_DEPOT_PATH/compiled/v$JULIA_VERSION"
if [ -d "$CACHE_DIR" ]; then
    rm -rf "$CACHE_DIR"
    echo "✓ Cache cleared: $CACHE_DIR"
else
    echo "⊘ No cache to clear"
fi
echo ""

# Step 1: Remove visualization packages (headless incompatible)
printf '=%.0s' {1..70}; echo ""
echo "STEP 1: Removing visualization packages..."
printf '=%.0s' {1..70}; echo ""
julia --project=. scripts/setup_headless.jl

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to remove visualization packages!"
    exit 1
fi

# Step 2: Verify visualization packages were removed
printf '=%.0s' {1..70}; echo ""
echo "STEP 2: Verifying headless setup..."
printf '=%.0s' {1..70}; echo ""
julia --project=. -e '
using Pkg
deps = keys(Pkg.project().dependencies)
bad_packages = ["GLMakie", "FileIO", "ColorSchemes", "LaTeXStrings", "MAT", "HDF5"]
found = filter(p -> p in deps, bad_packages)
if !isempty(found)
    println(stderr, "ERROR: Visualization packages still present in Project.toml!")
    println(stderr, "Found: ", join(found, ", "))
    println(stderr, "This will cause precompilation issues on headless HPC.")
    println(stderr, "Run scripts/setup_headless.jl manually to fix.")
    exit(1)
end
println("✓ No headless-incompatible packages found")
'

if [ $? -ne 0 ]; then
    exit 1
fi

# Step 3: Install remaining packages (including JLD2)
printf '=%.0s' {1..70}; echo ""
echo "STEP 3: Installing core packages..."
printf '=%.0s' {1..70}; echo ""
julia --project=. -e 'using Pkg; Pkg.instantiate()'

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install packages!"
    exit 1
fi

# Step 4: Configure MPI
printf '=%.0s' {1..70}; echo ""
echo "STEP 4: Configuring MPI..."
printf '=%.0s' {1..70}; echo ""
julia --project=. scripts/setup_mpi.jl

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to configure MPI!"
    exit 1
fi

# Step 5: Test
printf '=%.0s' {1..70}; echo ""
echo "STEP 5: Testing HyQMOM..."
printf '=%.0s' {1..70}; echo ""
julia --project=. -e 'using HyQMOM; println("✓ HyQMOM loaded successfully!")'

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to load HyQMOM!"
    exit 1
fi

# Summary
echo ""
printf '=%.0s' {1..70}; echo ""
echo "✓ SETUP COMPLETE!"
printf '=%.0s' {1..70}; echo ""
echo "Julia depot: $JULIA_DEPOT_PATH"
echo "GLMakie: REMOVED (headless)"
echo "JLD2: INSTALLED (snapshots)"
echo "MPI: CONFIGURED (system)"
echo ""
echo "Next: Submit simulation jobs using slurm/hyqmom_base.sbatch"
printf '=%.0s' {1..70}; echo ""