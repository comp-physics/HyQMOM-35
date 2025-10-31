# HyQMOM Quick Start for HPC

Ultra-condensed guide for running HyQMOM on headless HPC clusters.

## One-Time Setup (5 minutes)

```bash
# 1. Clone and checkout
git clone https://github.com/comp-physics/HyQMOM-35.git
cd HyQMOM-35
git checkout better-hpc
cd HyQMOM.jl

# 2. Load modules
module load julia/1.11
module load openmpi/4.1.5  # or your MPI version

# 3. Set Julia depot to exec-capable filesystem
export JULIA_DEPOT_PATH="$HOME/scratch/.julia"  # or /scratch/$USER/.julia
mkdir -p $JULIA_DEPOT_PATH

# 4. Remove headless-incompatible packages
julia --project=. scripts/setup_headless.jl

# 5. Install packages
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 6. Configure MPI
julia --project=. scripts/setup_mpi.jl

# 7. Test
julia --project=. -e 'using HyQMOM; println("✓ Ready!")'
```

**What this does:**
- ✅ Removes GLMakie, FileIO, ColorSchemes, LaTeXStrings (OpenGL/X11)
- ✅ Removes MAT, HDF5 (cause MPI JLL conflicts)
- ✅ Keeps JLD2 (snapshot I/O)
- ✅ Adds MPIPreferences and configures system MPI

**Result:** Your `Project.toml` will have only core deps (no visualization).

## Running Simulations

### Quick Test Run

```bash
# Interactive session or script
module load julia/1.11 openmpi/4.1.5
export JULIA_DEPOT_PATH="$HOME/scratch/.julia"
export HYQMOM_SKIP_PLOTTING=true
export OMPI_MCA_btl=^openib  # Silence InfiniBand warnings

cd HyQMOM.jl
srun --nodes=1 --ntasks=8 --mpi=pmix \
  julia --project=. examples/run_3d_custom_jets.jl \
  --Nx 40 --Ny 40 --Nz 40 --tmax 0.01 --snapshot-interval 10 --no-viz
```

### Production Run (Slurm)

Edit `scripts/example_slurm_run.sh`:
- Update `#SBATCH` directives for your cluster
- Set simulation parameters (grid size, Ma, Kn, etc.)

Submit:
```bash
sbatch scripts/example_slurm_run.sh
```

## Minimal Environment Variables

```bash
# Required
export JULIA_DEPOT_PATH="$HOME/scratch/.julia"  # exec-capable path
export HYQMOM_SKIP_PLOTTING=true                # skip viz at runtime

# Recommended (avoid warnings)
export OMPI_MCA_btl=^openib                     # disable openib if not configured
export JULIA_PKG_PRECOMPILE_AUTO=0              # no auto-precompile in jobs
export JULIA_NUM_THREADS=1                      # single-threaded Julia
export OMP_NUM_THREADS=1                        # no OpenMP
```

## Viewing Results

### Copy snapshots to local machine

```bash
# On local machine
scp user@cluster:/path/to/HyQMOM.jl/snapshots_*.jld2 .

# Visualize (GLMakie works on desktop)
cd /path/to/local/HyQMOM.jl
julia --project=. visualize_jld2.jl snapshots_*.jld2
```

## Troubleshooting

### "GLMakie failed to precompile"
You ran `Pkg.instantiate()` before `setup_headless.jl`. Fix:
```bash
julia --project=. scripts/setup_headless.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### "undefined symbol: ompi_instance_count"
MAT/HDF5 still present (pulls OpenMPI_jll). Fix:
```bash
julia --project=. -e 'using Pkg; Pkg.rm(["MAT","HDF5"]); Pkg.resolve()'
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### "LocalPreferences.toml not created"
MPI module not loaded when running setup_mpi.jl. Fix:
```bash
module load openmpi/4.1.5
julia --project=. scripts/setup_mpi.jl
```

### InfiniBand/UCX warnings
Add to your run script before `srun`:
```bash
export OMPI_MCA_btl=^openib
export UCX_WARN_UNUSED_ENV_VARS=n  # optional: suppress UCX warnings
```

### Jobs fail silently
Check your Julia depot is on an **executable** filesystem:
```bash
mount | grep $(dirname $JULIA_DEPOT_PATH)
# Should NOT show "noexec"
```

## What's Different from Desktop?

| Feature | Desktop | HPC Headless |
|---------|---------|--------------|
| **Visualization** | ✅ GLMakie | ❌ Removed |
| **PNG export** | ✅ FileIO | ❌ Removed |
| **MAT files** | ✅ MAT.jl | ❌ Removed (MPI conflict) |
| **JLD2 snapshots** | ✅ Works | ✅ **Works** |
| **MPI parallel** | ✅ Works | ✅ **Works** |
| **Interactive viz** | ✅ Live | ❌ Copy files to desktop |

## Command-Line Arguments

```bash
--Nx, --Ny, --Nz <N>     # Grid resolution (default: 40, 40, 20)
--tmax <time>            # Simulation end time (default: 0.05)
--Ma <mach>              # Mach number (default: 0.0)
--Kn <knudsen>           # Knudsen number (default: 1.0)
--CFL <number>           # CFL number (default: 0.7)
--snapshot-interval <N>  # Save every N steps (0=disabled)
--no-viz                 # Disable visualization (flag, no value)
--config <name>          # Jet config: crossing, triple-jet, quad-jet, etc.
```

Example:
```bash
julia --project=. examples/run_3d_custom_jets.jl \
  --Nx 80 --Ny 80 --Nz 80 --Ma 70.0 --Kn 10.0 \
  --tmax 0.0025 --snapshot-interval 50 --no-viz
```

## Files Modified by Setup

After running `setup_headless.jl`, your **local** `Project.toml` will show:

**Removed:**
- GLMakie, FileIO, ColorSchemes, LaTeXStrings (visualization)
- MAT, HDF5 (MPI JLL conflicts)

**Kept:**
- JLD2 (snapshot I/O)
- MPI (parallelization)
- LinearAlgebra, StaticArrays, Printf (core)

**Added by `setup_mpi.jl`:**
- MPIPreferences (MPI configuration)

⚠️ **Don't commit these changes!** The modified `Project.toml` is for your HPC cluster only. Desktop/CI needs the full deps.

## Performance Tips

- Use `--snapshot-interval` ≥ 50 for large grids (avoid I/O overhead)
- Set `--CFL 0.5` for stability with high Ma
- Monitor: `squeue -u $USER` and `tail -f slurm-<jobid>.out`
- Typical runtime: 80³ grid, Ma=70, Kn=10, tmax=0.0025 → ~1-2 hours on 8 nodes

## Support

- Full guide: `HyQMOM.jl/README_HPC.md` (if exists, else see repo README)
- Issues: https://github.com/comp-physics/HyQMOM-35/issues
- Example scripts: `HyQMOM.jl/scripts/example_slurm_*.sh`

