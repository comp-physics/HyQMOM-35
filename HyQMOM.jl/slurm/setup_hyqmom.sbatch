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

# Step 1: Remove visualization packages (pure shell, no Julia at all)
echo ""
printf '=%.0s' {1..70}; echo ""
echo "STEP 1: Removing visualization packages..."
printf '=%.0s' {1..70}; echo ""

# Backup Project.toml
cp -f Project.toml Project.toml.bak 2>/dev/null || true

# Packages to remove from [deps] and [compat] sections
PKGS=("GLMakie" "FileIO" "ColorSchemes" "LaTeXStrings" "MAT" "HDF5")

for P in "${PKGS[@]}"; do
  # Remove from [deps] (handles with or without quotes)
  if awk '/^\[deps\]/{flag=1;next} /^\[/{flag=0} flag && $0 ~ "^[[:space:]]*\"?'$P'\"?[[:space:]]*=" {found=1} END{exit !found}' Project.toml; then
    sed -i.tmp '/^\[deps\]/,/^\[/{/^[[:space:]]*"?'$P'"?[[:space:]]*=/d}' Project.toml
    rm -f Project.toml.tmp
    echo "  ✓ Removed from [deps]: $P"
  else
    echo "  ⊘ Skipped [deps]: $P (not found)"
  fi

  # Remove from [compat] (handles with or without quotes)
  if awk '/^\[compat\]/{flag=1;next} /^\[/{flag=0} flag && $0 ~ "^[[:space:]]*\"?'$P'\"?[[:space:]]*=" {found=1} END{exit !found}' Project.toml; then
    sed -i.tmp '/^\[compat\]/,/^\[/{/^[[:space:]]*"?'$P'"?[[:space:]]*=/d}' Project.toml
    rm -f Project.toml.tmp
    echo "  ✓ Removed from [compat]: $P"
  fi
done

# Debug: show any lingering compat entries for these packages
echo ""
echo "Checking for lingering compat entries..."
LINGERING=$(awk '/^\[compat\]/{flag=1;next} /^\[/{flag=0} flag {print}' Project.toml | grep -E '^[[:space:]]*(\"?(GLMakie|FileIO|ColorSchemes|LaTeXStrings|MAT|HDF5)\"?[[:space:]]*=)' || true)
if [ -n "$LINGERING" ]; then
  echo "⚠ WARNING: Found lingering compat entries:"
  echo "$LINGERING"
else
  echo "  ✓ No lingering entries"
fi

# Verify nothing remains in [deps]
if awk '/^\[deps\]/{flag=1;next} /^\[/{flag=0} flag && $0 ~ "^(GLMakie|FileIO|ColorSchemes|LaTeXStrings|MAT|HDF5)[[:space:]]*=" {exit 1}' Project.toml; then
  echo ""
  echo "✓ Project.toml updated"
else
  echo ""
  echo "ERROR: Some headless-incompatible packages still present in [deps]" >&2
  echo "Hint: compare Project.toml and Project.toml.bak" >&2
  exit 1
fi

# Clear entire compiled cache after package removal
echo ""
echo "Clearing compiled cache..."
if [ -d "$JULIA_DEPOT_PATH/compiled" ]; then
  rm -rf "$JULIA_DEPOT_PATH/compiled"
  echo "✓ Entire cache cleared"
else
  echo "⊘ No cache to clear"
fi
echo ""

# Step 2: Verify visualization packages were removed (pure shell, no Julia)
printf '=%.0s' {1..70}; echo ""
echo "STEP 2: Verifying headless setup..."
printf '=%.0s' {1..70}; echo ""

# Check Project.toml [deps] section for bad packages
if awk '/^\[deps\]/{flag=1;next} /^\[/{flag=0} flag && $0 ~ "^(GLMakie|FileIO|ColorSchemes|LaTeXStrings|MAT|HDF5)[[:space:]]*=" {exit 1}' Project.toml; then
  echo "✓ No headless-incompatible packages found in Project.toml"
  echo ""
else
  echo "ERROR: Headless-incompatible packages still in Project.toml [deps]!" >&2
  echo "Check Project.toml manually and compare with Project.toml.bak" >&2
  echo ""
  exit 1
fi

# Step 3: Clean slate - wipe Julia depot for fresh rebuild
echo "Wiping Julia depot for clean rebuild..."
if [ -d "$JULIA_DEPOT_PATH" ]; then
  # Backup registries if they exist (saves re-download time)
  if [ -d "$JULIA_DEPOT_PATH/registries" ]; then
    mv "$JULIA_DEPOT_PATH/registries" "$JULIA_DEPOT_PATH/../registries.bak" 2>/dev/null || true
  fi
  
  # Remove entire depot (nuclear option for depot corruption)
  rm -rf "$JULIA_DEPOT_PATH"
  mkdir -p "$JULIA_DEPOT_PATH"
  
  # Restore registries
  if [ -d "$JULIA_DEPOT_PATH/../registries.bak" ]; then
    mv "$JULIA_DEPOT_PATH/../registries.bak" "$JULIA_DEPOT_PATH/registries"
    echo "✓ Julia depot wiped (registries preserved)"
  else
    echo "✓ Julia depot wiped (will re-download registries)"
  fi
else
  echo "⊘ No depot to wipe"
fi
echo ""

# Delete Manifest.toml for clean rebuild
if [ -f "Manifest.toml" ]; then
  rm -f Manifest.toml
  echo "✓ Removed Manifest.toml"
fi
echo ""

# Step 4: Install remaining packages (including JLD2)
printf '=%.0s' {1..70}; echo ""
echo "STEP 4: Installing core packages..."
printf '=%.0s' {1..70}; echo ""
julia --project=. -e 'using Pkg; Pkg.instantiate()'

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install packages!"
    exit 1
fi

# Step 5: Configure MPI
printf '=%.0s' {1..70}; echo ""
echo "STEP 5: Configuring MPI..."
printf '=%.0s' {1..70}; echo ""
julia --project=. scripts/setup_mpi.jl

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to configure MPI!"
    exit 1
fi

# Step 6: Test
printf '=%.0s' {1..70}; echo ""
echo "STEP 6: Testing HyQMOM..."
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