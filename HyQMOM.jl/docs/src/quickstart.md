# Quickstart Guide

Get up and running with HyQMOM.jl in just a few minutes!

## Installation

First, clone the repository and set up the Julia environment:

```bash
git clone https://github.com/comp-physics/HyQMOM.jl.git
cd HyQMOM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This will download and install all required dependencies.

## Your First Simulation

Run an interactive 3D visualization example:

```bash
julia --project=. examples/run_3d_jets_timeseries.jl
```

This will:
1. Run a 3D crossing jets simulation with default parameters
2. Launch an interactive GLMakie viewer
3. Show real-time 3D isosurfaces of density and velocity fields
4. Save snapshots to a `.jld2` file for later analysis

### What You'll See

The interactive viewer includes:
- **Time slider**: Step through simulation snapshots
- **Play/Pause/Reset**: Animate the time evolution
- **Quantity buttons**: Switch between density, U/V/W velocities, pressure, temperature
- **Isosurface sliders**: Adjust visualization levels
- **Mouse controls**: Drag to rotate, scroll to zoom

## Quick Parameter Adjustments

Try different simulation parameters:

```bash
# Quick low-resolution test (faster)
julia --project=. examples/run_3d_jets_timeseries.jl --Np 20 --tmax 0.01

# Higher resolution (more detailed)
julia --project=. examples/run_3d_jets_timeseries.jl --Np 60 --tmax 0.1

# Different physics parameters
julia --project=. examples/run_3d_jets_timeseries.jl --Ma 1.5 --Kn 0.5
```

### Key Parameters
- `--Np N`: Grid resolution in x,y directions (default: 40)
- `--Nz N`: Grid resolution in z direction (default: 40)  
- `--tmax T`: Maximum simulation time (default: 0.2)
- `--Ma M`: Mach number (default: 1.0)
- `--Kn K`: Knudsen number (default: 1.0)
- `--CFL C`: CFL number for stability (default: 0.7)

Use `--help` to see all available options.

## MPI Parallel Execution

Scale up with MPI for larger simulations:

```bash
# Run with 4 MPI processes
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100

# High-resolution production run
mpiexec -n 8 julia --project=. examples/run_3d_jets_timeseries.jl --Np 120 --Nz 60
```

The visualization will automatically appear on rank 0, while all ranks participate in the computation.

## Visualizing Saved Results

After running simulations, you can re-visualize the saved `.jld2` files:

```bash
# Auto-detect and visualize .jld2 files in current directory
julia visualize_jld2.jl

# Or specify a file directly
julia visualize_jld2.jl snapshots_crossing_Ma1.0_t0.3_N30.jld2
```

## Headless/HPC Systems

For systems without display capabilities (compute clusters, remote servers):

```bash
# Skip visualization but save snapshots
julia --project=. examples/run_3d_jets_timeseries.jl --no-viz true

# Or set environment variable
export HYQMOM_SKIP_PLOTTING=true
julia --project=. examples/run_3d_jets_timeseries.jl
```

Transfer the `.jld2` files to a system with display capabilities and use `julia visualize_jld2.jl` to view them.

## Next Steps

- **User Guide**: Comprehensive parameter documentation and workflows
- **MPI & Parallelization**: Detailed parallel computing guide
- **Tutorials**: Step-by-step examples
- **API Reference**: Complete function documentation

## Troubleshooting

### GLMakie doesn't open
```bash
# Test GLMakie installation
julia --project=. -e 'using GLMakie'

# On remote systems, enable X11 forwarding
ssh -Y user@host
```

### Out of memory
```bash
# Reduce resolution
julia ... --Np 30 --Nz 15

# Increase snapshot interval or disable
julia ... --snapshot-interval 10
julia ... --snapshot-interval 0
```

### Simulation crashes (NaN values)
```bash
# Reduce CFL number for stability
julia ... --CFL 0.5

# Reduce Mach number
julia ... --Ma 0.7

# Increase resolution
julia ... --Np 60
```
