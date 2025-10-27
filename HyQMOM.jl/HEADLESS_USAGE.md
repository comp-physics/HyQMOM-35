# Running HyQMOM on Headless Systems (No Display/X11)

This guide explains how to run simulations on compute clusters, HPC systems, or remote servers without display capabilities, then visualize results later on your local machine.

## Quick Start

### 1. Run Simulation on Headless System

Use the `--no-viz` flag to disable visualization while still saving `.jld2` snapshot files:

```bash
# Serial run
julia --project=. examples/run_3d_custom_jets.jl \
  --config crossing \
  --snapshot-interval 5 \
  --no-viz true

# MPI parallel run (recommended for HPC)
mpiexec -n 16 julia --project=. examples/run_3d_custom_jets.jl \
  --config quad-jet \
  --Nx 80 --Ny 80 --Nz 40 \
  --tmax 0.2 \
  --snapshot-interval 10 \
  --no-viz true
```

### 2. Transfer Results to Local Machine

```bash
# From local machine
scp user@cluster:~/path/to/snapshots_*.jld2 .
```

### 3. Visualize Locally

```bash
# Interactive 3D viewer
julia visualize_jld2.jl snapshots_crossing_Ma1.0_t0.2_N80.jld2

# Or let it auto-detect
julia visualize_jld2.jl
```

## What Gets Saved

With `--no-viz true`, the simulation will:
- Run normally (no performance impact)
- Save `.jld2` files with all snapshots
- Include standardized moments (S) automatically
- Work with both serial and MPI runs
- Skip launching the interactive viewer
- Not load GLMakie at all (faster startup, less memory)

## Complete Workflow Examples

### Example 1: Quick Test Run

```bash
# On HPC cluster
julia --project=. examples/run_3d_custom_jets.jl \
  --config triple-jet \
  --Nx 30 --Ny 30 --Nz 20 \
  --tmax 0.1 \
  --snapshot-interval 5 \
  --no-viz true

# Output: snapshots_triple-jet_Ma0.0_t0.1_N30.jld2

# On local machine
scp cluster:snapshots_triple-jet_Ma0.0_t0.1_N30.jld2 .
julia visualize_jld2.jl
```

### Example 2: Production MPI Run

```bash
# On HPC cluster (submit as batch job)
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00

cd $WORK/HyQMOM.jl
module load julia/1.10
module load openmpi

mpiexec -n 32 julia --project=. examples/run_3d_custom_jets.jl \
  --config crossing \
  --Nx 120 --Ny 120 --Nz 60 \
  --tmax 0.5 \
  --Ma 1.2 \
  --Kn 0.8 \
  --CFL 0.6 \
  --snapshot-interval 20 \
  --no-viz true

# Output: snapshots_crossing_Ma1.2_t0.5_N120.jld2
```

### Example 3: Parameter Sweep

```bash
# On HPC: run_sweep.sh
#!/bin/bash

for Ma in 0.5 1.0 1.5 2.0; do
  for Nx in 40 60 80; do
    echo "Running Ma=$Ma, Nx=$Nx"
    
    julia --project=. examples/run_3d_custom_jets.jl \
      --config crossing \
      --Nx $Nx --Ny $Nx --Nz $((Nx/2)) \
      --Ma $Ma \
      --tmax 0.2 \
      --snapshot-interval 10 \
      --no-viz true
      
    echo "Saved: snapshots_crossing_Ma${Ma}_t0.2_N${Nx}.jld2"
  done
done

# Download all results
# rsync -avz cluster:~/HyQMOM.jl/snapshots_*.jld2 ./results/
```

## Snapshot File Contents

Each `.jld2` file contains:
- `snapshots`: Vector of snapshots with fields:
  - `M`: Raw moment field (Nx x Ny x Nz x 35)
  - `t`: Simulation time
  - `step`: Time step number
  - `S`: Standardized moments (always saved)
  - `C`: Central moments (always saved)
- `grid`: Grid coordinates (xm, ym, zm, dx, dy, dz)
- `params`: Simulation parameters
- `params_with_ic`: Parameters including initial conditions

## Visualization Options

### Interactive 3D Time-Series (Recommended)
```julia
using HyQMOM, JLD2, GLMakie
# The viewer shows physical space (left) and moment space (middle panel)
# Moment space displays (S110, S101, S011) if S field is saved
interactive_3d_timeseries_streaming("snapshots_file.jld2", grid, params_with_ic)
```

### Static PyPlot Figures (for publications)
```julia
using PyPlot
M = snapshots[end].M
Nx = params.Nx; Nz = params.Nz
plot_final_results(M, grid.xm, grid.ym, Nx, 35; 
                   zm=grid.zm, Nz=Nz, save_figures=true)
```

## Troubleshooting

### Error: "cannot open display"
**Solution:** Use `--no-viz true` flag


### Cannot visualize on local machine
**Cause:** GLMakie not installed or display issues

**Solution:**
```bash
cd HyQMOM.jl
julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'
julia visualize_jld2.jl
```


## Environment Variables

For package-level control (affects all scripts):

```bash
# Skip loading GLMakie entirely (for CI/testing)
export HYQMOM_SKIP_PLOTTING=true

# Or use CI flag (auto-detected)
export CI=true
```

These environment variables prevent GLMakie from loading at the package level, making the entire package lightweight. Use these for:
- Continuous Integration
- Unit testing
- When GLMakie cannot be installed

For running simulations and saving results, prefer the `--no-viz` flag instead.

