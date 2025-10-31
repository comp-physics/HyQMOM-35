# HyQMOM Setup for HPC Clusters

This guide provides instructions for setting up and running HyQMOM on headless HPC systems (no X11/OpenGL).

## Quick Setup for Headless Systems

### 1. Clone and Checkout Branch

```bash
git clone https://github.com/comp-physics/HyQMOM-35.git
cd HyQMOM-35
git checkout better-hpc
cd HyQMOM.jl
```

### 2. One-Time Setup Script

Create `setup_hyqmom.sh`:

```bash
#!/bin/bash
#SBATCH -J hyqmom_setup
#SBATCH -A your-account
#SBATCH -p partition-name
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=16GB

# Load modules
module load julia/1.11
module load openmpi/4.1.5  # or your MPI version

# Use exec-capable depot (CRITICAL - not HOME if it's noexec!)
export JULIA_DEPOT_PATH="/scratch/$USER/.julia"
mkdir -p $JULIA_DEPOT_PATH

cd HyQMOM.jl

# Step 1: Remove visualization packages (headless incompatible)
echo "Step 1: Removing visualization packages..."
julia --project=. scripts/setup_headless.jl

# Step 2: Install remaining packages (including JLD2)
echo ""
echo "Step 2: Installing core packages..."
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Step 3: Configure MPI
echo ""
echo "Step 3: Configuring MPI..."
julia --project=. scripts/setup_mpi.jl

# Step 4: Test
echo ""
echo "Step 4: Testing..."
julia --project=. -e 'using HyQMOM; println("✓ HyQMOM loaded successfully!")'

echo ""
echo "✓ HyQMOM setup complete!"
```

Run setup:
```bash
sbatch setup_hyqmom.sh
```

### 3. Manual Setup (Alternative)

If you prefer to run commands interactively:

```bash
# Load modules
module load julia/1.11
module load openmpi/4.1.5

# Set Julia depot
export JULIA_DEPOT_PATH="/scratch/$USER/.julia"
mkdir -p $JULIA_DEPOT_PATH

cd HyQMOM.jl

# Remove visualization packages
julia --project=. scripts/setup_headless.jl

# Install packages
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Configure MPI
julia --project=. scripts/setup_mpi.jl

# Test
julia --project=. -e 'using HyQMOM; println("✓ Ready!")'
```

## Running Simulations

### Production Run Script

Create `run_simulation.sh`:

```bash
#!/bin/bash
#SBATCH -J hyqmom_Ma70_Kn10
#SBATCH -A your-account
#SBATCH -p partition-name
#SBATCH -N 8
#SBATCH --ntasks-per-node=96
#SBATCH --time=8:00:00
#SBATCH --mem=128GB

# Load modules (SAME as setup)
module load julia/1.11
module load openmpi/4.1.5

# Use SAME depot as setup
export JULIA_DEPOT_PATH="/scratch/$USER/.julia"

# Disable plotting and precompilation
export HYQMOM_SKIP_PLOTTING=true
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OMP_NUM_THREADS=1

cd HyQMOM.jl

# Run simulation with MPI
srun --mpi=pmix julia --project=. examples/run_3d_custom_jets.jl \
  --Nx 80 --Ny 80 --Nz 80 \
  --tmax 0.0025 \
  --Ma 70.0 \
  --Kn 10.0 \
  --snapshot-interval 50 \
  --no-viz true

echo ""
echo "✓ Simulation complete!"
echo "  Output files: snapshots_*.jld2"
```

Submit:
```bash
sbatch run_simulation.sh
```

## What Gets Installed

### ✅ Installed (Always)
- **JLD2** - Snapshot file I/O (`.jld2` format)
- **MPI** - Parallel execution
- **LinearAlgebra, StaticArrays, Printf** - Core numerics
- **Dates, DelimitedFiles, MAT** - Utilities

### ❌ Removed (Headless Incompatible)
- **GLMakie** - Interactive 3D visualization (requires OpenGL)
- **FileIO** - PNG export from viewer (visualization only)
- **ColorSchemes** - Color maps for plots
- **LaTeXStrings** - LaTeX labels in plots

## Viewing Results

### Copy Snapshots to Local Machine

After job completes:

```bash
# On your local machine
scp username@cluster:/path/to/HyQMOM.jl/snapshots_*.jld2 .

# View locally (where GLMakie IS available)
cd /path/to/local/HyQMOM.jl
julia --project=. visualize_jld2.jl snapshots_*.jld2
```

## Troubleshooting

### "GLMakie failed to precompile"

This is expected! The setup scripts remove GLMakie before installation. If you see this error, it means you ran `Pkg.instantiate()` **before** running `setup_headless.jl`.

**Fix:**
```bash
cd HyQMOM.jl
julia --project=. scripts/setup_headless.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### "MPIPreferences not found"

Run the MPI setup script:
```bash
julia --project=. scripts/setup_mpi.jl
```

### "Permission denied" for libmpi.so

Your Julia depot is on a `noexec` filesystem. Set it to an exec-capable location:
```bash
export JULIA_DEPOT_PATH="/scratch/$USER/.julia"  # or /tmp/$USER/.julia
```

### LocalPreferences.toml not created

Check that you:
1. Loaded the MPI module (`module load openmpi/...`)
2. Ran `setup_mpi.jl` script
3. Check file exists: `cat HyQMOM.jl/LocalPreferences.toml`

Should contain:
```toml
[MPIPreferences]
binary = "system"
```

### "Being precompiled by another process"

This means multiple MPI ranks are trying to precompile simultaneously. This shouldn't happen if you:
1. Ran setup on a single node first
2. Set `JULIA_PKG_PRECOMPILE_AUTO=0` in run script

If it persists, clear depot and re-run setup:
```bash
rm -rf /scratch/$USER/.julia
# Re-run setup
```

## Performance Notes

### Expected Startup Times

- **After Setup:** < 5 seconds to load HyQMOM
- **Before Setup (with GLMakie):** 5-10 minutes per job

### Snapshot File Sizes

Typical sizes for `--snapshot-interval 50`:
- **80x80x80 grid:** ~50-100 MB per file
- **160x160x160 grid:** ~400-800 MB per file

Snapshots include:
- Moment field `M` (35 components)
- Standardized moments `S` (4 components: S110, S101, S011)
- Central moments `C` (4 components)
- Time and step number

## Environment Variables Reference

| Variable | Purpose | Recommended Value |
|----------|---------|-------------------|
| `JULIA_DEPOT_PATH` | Package storage | `/scratch/$USER/.julia` |
| `HYQMOM_SKIP_PLOTTING` | Skip viz packages | `true` |
| `JULIA_PKG_PRECOMPILE_AUTO` | Disable auto-precompile | `0` |
| `JULIA_NUM_THREADS` | Julia threads | `1` |
| `OMP_NUM_THREADS` | OpenMP threads | `1` |

## Command-Line Arguments

### run_3d_custom_jets.jl

```bash
--Nx <N>              # Grid points in x (default: 20)
--Ny <N>              # Grid points in y (default: 20)
--Nz <N>              # Grid points in z (default: 20)
--tmax <time>         # Simulation end time (default: 0.1)
--Ma <mach>           # Mach number (default: 0.0)
--Kn <knudsen>        # Knudsen number (default: 1.0)
--CFL <number>        # CFL number (default: 0.3)
--snapshot-interval <N>  # Save every N steps (default: 1, 0=disable)
--no-viz <bool>       # Disable visualization (default: false)
--config <name>       # Jet configuration: crossing, crossing2D, triple-jet, quad-jet, etc.
```

### Examples

```bash
# Small test run
srun julia --project=. examples/run_3d_custom_jets.jl \
  --Nx 20 --Ny 20 --Nz 20 --tmax 0.01 --snapshot-interval 10 --no-viz true

# Production run
srun julia --project=. examples/run_3d_custom_jets.jl \
  --Nx 160 --Ny 160 --Nz 160 --tmax 0.005 --Ma 70.0 --Kn 10.0 \
  --snapshot-interval 100 --no-viz true

# Different jet configuration
srun julia --project=. examples/run_3d_custom_jets.jl \
  --config quad-jet --Nx 80 --Ny 80 --Nz 80 --no-viz true
```

## Support

For issues specific to HPC setup, check:
1. This README
2. GitHub Issues: https://github.com/comp-physics/HyQMOM-35/issues
3. Main README: `../README.md`

For general HyQMOM usage, see:
- `../README.md` - Main documentation
- `examples/` - Example scripts
- `docs/` - Full documentation (if available)

