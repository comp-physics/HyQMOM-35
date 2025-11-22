# HyQMOM.jl

[![CI](https://github.com/comp-physics/HyQMOM.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/comp-physics/HyQMOM.jl/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/hyqmomjl/badge/?version=latest)](https://hyqmomjl.readthedocs.io/en/latest/)
[![codecov](https://codecov.io/gh/comp-physics/HyQMOM.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/comp-physics/HyQMOM.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17682195.svg)](https://doi.org/10.5281/zenodo.17682195)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia Version](https://img.shields.io/badge/julia-v1.9+-blue.svg)](https://julialang.org/downloads/)

3D Hyperbolic Quadrature Method of Moments (HyQMOM) solver for the Boltzmann equation with BGK collision operator, featuring MPI parallelization and interactive visualization.

## Quick Start

### Installation

```bash
# Clone repository (if not already done)
git clone https://github.com/comp-physics/HyQMOM.jl.git
cd HyQMOM.jl

# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run Your First Simulation

```bash
# Interactive 3D visualization (recommended)
julia --project=. examples/run_3d_jets_timeseries.jl

# With MPI parallelization
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 60
```

## Features

- **3D moment-based kinetic solver** for Boltzmann-BGK equation
- **MPI parallelization** with domain decomposition
- **Interactive 3D visualization** with GLMakie (time-series animation)
- **Flexible parameter control** via command-line arguments
- **Serial or parallel** execution (automatic detection)
- **Snapshot collection** for time evolution analysis

## Usage

### Basic Examples

```bash
# Quick low-resolution test
julia --project=. examples/run_3d_jets_timeseries.jl --Np 20 --tmax 0.01

# Production run with 4 MPI ranks
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100 --tmax 0.5

# Custom domain and physics
julia --project=. examples/run_3d_jets_timeseries.jl \
  --Np 60 --tmax 0.1 --Ma 1.5 --Kn 0.5 --CFL 0.6
```

### Common Parameters

```bash
--Np N             # Grid resolution x,y (default: 40)
--Nz N             # Grid resolution z (default: 40)
--tmax T           # Maximum simulation time (default: 0.2)
--Ma M             # Mach number (default: 1.0)
--Kn K             # Knudsen number (default: 1.0)
--CFL C            # CFL number (default: 0.7)
--snapshot-interval N  # Save every N steps for visualization
--help             # Show all available options
```

See `examples/README.md` for complete parameter documentation.

## Examples

### Visualization 

**`examples/run_3d_jets_timeseries.jl`** - Interactive 3D time-series viewer

```bash
julia --project=. examples/run_3d_jets_timeseries.jl
```

Features:
- Real-time 3D isosurface visualization
- Time slider with play/pause animation  
- Multiple quantities (density, U/V/W velocities, pressure, temperature)
- Works in both serial and MPI parallel modes

### Custom Initial Conditions 

**`examples/run_3d_custom_jets.jl`** - Flexible jet configurations

```bash
# Triple jet configuration
julia --project=. examples/run_3d_custom_jets.jl --config triple-jet

# Four jets converging to center
julia --project=. examples/run_3d_custom_jets.jl --config quad-jet

# Create your own configurations!
```

Features:
- Predefined configurations: crossing, triple-jet, quad-jet, vertical-jet, spiral
- Fully customizable: specify center, width, velocity for each cubic region
- Easy to add new configurations

See `examples/CUSTOM_INITIAL_CONDITIONS.md` for detailed documentation.

### Custom Configurations

**`examples/run_3d_custom_jets.jl`** - Flexible jet configurations

```bash
julia --project=. examples/run_3d_custom_jets.jl --config crossing
julia --project=. examples/run_3d_custom_jets.jl --config triple-jet
```

Choose from multiple predefined configurations or create your own.

## Visualization

### Quick Visualization of Results

After running simulations with `examples/run_3d_custom_jets.jl`, you'll get `.jld2` snapshot files. Use the provided script to quickly visualize them:

```bash
# Auto-detect and visualize .jld2 files
julia visualize_jld2.jl

# Or specify a file directly
julia visualize_jld2.jl snapshots_crossing_Ma1.0_t0.3_N30.jld2
```

The script will:
- **Auto-install** GLMakie if needed
- **Auto-find** `.jld2` files in the current directory
- **Launch** interactive 3D time-series viewer
- **Show** both physical space (density/velocity isosurfaces) and moment space (if available)

### Interactive 3D Visualization

HyQMOM includes interactive 3D visualization using GLMakie. GLMakie is included as a dependency and works out of the box.

**Controls:**
- Time slider: Step through simulation snapshots
- Play/Pause/Reset: Animate the time evolution
- Quantity buttons: Switch between density, velocities, etc.
- Isosurface sliders: Adjust visualization levels
- Mouse: Drag to rotate, scroll to zoom

**Manual Usage:**
```julia
using HyQMOM, JLD2, GLMakie
@load "snapshots_file.jld2" snapshots grid params params_with_ic

# Interactive time-series viewer (streams from file)
# The middle panel shows standardized moment space (S110, S101, S011) if available
interactive_3d_timeseries_streaming("snapshots_file.jld2", grid, params_with_ic)
```

### Headless Systems (HPC/Clusters)

**For running on systems without display/X11** (compute clusters, HPC, remote servers):

```bash
# Use the --no-viz flag to skip visualization but still save .jld2 files
julia --project=. examples/run_3d_custom_jets.jl --no-viz true --snapshot-interval 5

# Or for run_3d_jets_timeseries.jl
julia --project=. examples/run_3d_jets_timeseries.jl --no-viz true

# View results later on a system with display:
julia visualize_jld2.jl
```

**See [HEADLESS_USAGE.md](HEADLESS_USAGE.md) for complete workflow examples, MPI batch jobs, and parameter sweeps.**

**For CI pipelines or when installing the package** (skip loading GLMakie entirely):

```bash
export HYQMOM_SKIP_PLOTTING=true
# or
export CI=true  # Automatically detected
```

This skips loading GLMakie (~100 packages), making the package lightweight for testing and headless systems.

**Example CI configuration (GitHub Actions):**

```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    env:
      HYQMOM_SKIP_PLOTTING: "true"
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v1
      - uses: julia-actions/cache@v1
      # Remove GLMakie to avoid X11/display errors in CI
      - run: julia --project=. -e 'using Pkg; Pkg.rm("GLMakie")'
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
```

See `.github/workflows/ci.yml` for the complete CI configuration.

## MPI Parallelization

All examples automatically support MPI parallelization:

```bash
# Serial (1 process)
julia --project=. examples/run_3d_jets_timeseries.jl --Np 40

# Parallel (4 processes)
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100

# Parallel (8 processes, high resolution)
mpiexec -n 8 julia --project=. examples/run_3d_jets_timeseries.jl --Np 120 --Nz 60
```

**Domain Decomposition:**
- x-y plane divided among MPI ranks
- z-direction: full data on all ranks (no decomposition)
- Automatic halo exchange for ghost cells
- Visualization only on rank 0

## API Usage

### Direct API

```julia
using HyQMOM
using MPI

MPI.Init()

# Set up parameters
params = (
    Np = 40,
    Nz = 40,
    tmax = 0.1,
    Ma = 0.0,
    Kn = 1.0,
    CFL = 0.7,
    # ... see examples/parse_params.jl for all options
)

# Run simulation with snapshot collection
snapshots, grid = run_simulation_with_snapshots(params; snapshot_interval=2)

# Or run without snapshots (just final result)
M_final, final_time, time_steps, grid = simulation_runner(params)

MPI.Finalize()
```

### With Visualization

```julia
using HyQMOM
import GLMakie

# Run simulation
snapshots, grid = run_simulation_with_snapshots(params; snapshot_interval=2)

# Launch interactive viewer (rank 0 only)
# Note: Save snapshots to file first, then use interactive_3d_timeseries_streaming
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @info "Save snapshots to file, then use: interactive_3d_timeseries_streaming(filename, grid, params)"
end
```

## Testing

```bash
# Run all tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run with MPI tests
cd test && bash run_mpi_tests.sh

# CI mode (skip visualization)
HYQMOM_SKIP_PLOTTING=true julia --project=. -e 'using Pkg; Pkg.test()'
```

## Performance Tips

### For Quick Tests
```bash
julia ... --Np 20 --Nz 10 --tmax 0.01 --snapshot-interval 10
```

### For Production
```bash
mpiexec -n 8 julia ... --Np 100 --Nz 50 --snapshot-interval 0
```

### Memory Optimization
- Reduce resolution: `--Np 30`
- Increase snapshot interval: `--snapshot-interval 10`
- Disable snapshots: `--snapshot-interval 0`

### Numerical Stability
- Reduce CFL: `--CFL 0.5`
- Reduce Mach number: `--Ma 0.7`
- Increase resolution: `--Np 60`

## Troubleshooting

### GLMakie viewer doesn't open
```bash
# Check GLMakie is installed
julia --project=. -e 'using GLMakie'

# On remote systems, enable X11 forwarding
ssh -Y user@host

# Or use custom jet configurations
julia --project=. examples/run_3d_custom_jets.jl --config crossing
```

### MPI errors
```bash
# Ensure MPI is properly configured
mpiexec --version

# Try fewer ranks
mpiexec -n 2 julia ...

# Check MPI.jl configuration
julia --project=. -e 'using MPI; println(MPI.versioninfo())'
```

### Out of memory
```bash
# Reduce resolution
julia ... --Np 30 --Nz 15

# Increase snapshot interval or disable
julia ... --snapshot-interval 10
julia ... --snapshot-interval 0

# Use more MPI ranks to distribute memory
mpiexec -n 8 julia ...
```

### NaN values / simulation crashes
```bash
# Reduce CFL number
julia ... --CFL 0.5

# Reduce Mach number
julia ... --Ma 0.7

# Increase resolution
julia ... --Np 60 --Nz 30
```

## Development

### Adding New Features

1. Core functionality goes in `src/`
2. Visualization in `src/visualization/`
3. Examples in `examples/`
4. Tests in `test/`

### Exported Functions

**Main simulation:**
- `simulation_runner(params)` - Run simulation, return final state
- `run_simulation_with_snapshots(params; snapshot_interval)` - Collect time snapshots

**Visualization:**
- `interactive_3d_timeseries_streaming(filename, grid, params)` - Streaming time-series viewer

**Core algorithms:**
- `hyqmom_3D`, `InitializeM4_35`, `Moments5_3D` - Moment operations
- `realizability`, `realizable_2D`, `realizable_3D` - Realizability checks
- `Flux_closure35_and_realizable_3D`, `flux_HLL`, `collision35` - Numerics

See `src/HyQMOM.jl` for complete list of exports.

## License

See `license.md` for licensing information.

---

**Quick Reference:**

```bash
# Basic usage
julia --project=. examples/run_3d_jets_timeseries.jl

# MPI parallel
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100

# Visualize results
julia visualize_jld2.jl

# Get help
julia --project=. examples/run_3d_jets_timeseries.jl --help

# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

