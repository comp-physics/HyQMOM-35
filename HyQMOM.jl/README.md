# HyQMOM.jl

3D Hyperbolic Quadrature Method of Moments (HyQMOM) solver for the Boltzmann equation with BGK collision operator, featuring MPI parallelization and interactive visualization.

## Quick Start

### Installation

```bash
cd HyQMOM.jl
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

### Interactive Visualization (Recommended)

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

### Static Plots

**`examples/run_3d_crossing_jets.jl`** - Static PyPlot visualization

```bash
julia --project=. examples/run_3d_crossing_jets.jl
```

Good for publications and batch processing.

## Visualization

### Interactive 3D Visualization

HyQMOM includes interactive 3D visualization using GLMakie. GLMakie is included as a dependency and works out of the box.

**Controls:**
- Time slider: Step through simulation snapshots
- Play/Pause/Reset: Animate the time evolution
- Quantity buttons: Switch between density, velocities, etc.
- Isosurface sliders: Adjust visualization levels
- Mouse: Drag to rotate, scroll to zoom

### Disabling Visualization (for CI/HPC)

For CI pipelines or HPC systems without graphics:

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
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    interactive_3d_timeseries(snapshots, grid, params)
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

## Project Structure

```
HyQMOM.jl/
├── src/
│   ├── HyQMOM.jl              # Main module
│   ├── moments/               # Moment operations (InitializeM4_35, hyqmom_3D, etc.)
│   ├── realizability/         # Realizability checks
│   ├── numerics/              # Flux closures, eigenvalues, collision operator
│   ├── mpi/                   # MPI utilities and halo exchange
│   ├── utils/                 # Utilities (moment indexing, diagnostics)
│   ├── visualization/         # GLMakie visualization (optional)
│   ├── simulation_runner.jl   # Main simulation loop
│   └── autogen/               # Auto-generated symbolic code
├── examples/
│   ├── run_3d_jets_timeseries.jl   # Interactive visualization example
│   ├── run_3d_crossing_jets.jl     # Static plots example
│   ├── parse_params.jl             # Parameter parsing utilities
│   └── README.md                    # Detailed examples documentation
├── test/                      # Test suite
├── Project.toml               # Package dependencies
└── README.md                  # This file
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
ssh -X user@host

# Or use static plots instead
julia --project=. examples/run_3d_crossing_jets.jl
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
- `interactive_3d_timeseries(snapshots, grid, params)` - Time-series viewer
- `interactive_3d_volume(M_final, grid, params)` - Single snapshot viewer

**Core algorithms:**
- `hyqmom_3D`, `InitializeM4_35`, `Moments5_3D` - Moment operations
- `realizability`, `realizable_2D`, `realizable_3D` - Realizability checks
- `Flux_closure35_and_realizable_3D`, `flux_HLL`, `collision35` - Numerics

See `src/HyQMOM.jl` for complete list of exports.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{hyqmom_jl,
  title = {HyQMOM.jl: 3D Hyperbolic Quadrature Method of Moments},
  author = {Spencer H. Bryngelson and contributors},
  year = {2024},
  url = {https://github.com/comp-physics/HyQMOM.jl}
}
```

## License

See `license.md` for licensing information.

## Contact

For questions, issues, or contributions:
- Open an issue on GitHub
- See `examples/README.md` for detailed usage examples
- Check the test suite in `test/` for API usage examples

---

**Quick Reference:**

```bash
# Basic usage
julia --project=. examples/run_3d_jets_timeseries.jl

# MPI parallel
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100

# Get help
julia --project=. examples/run_3d_jets_timeseries.jl --help

# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

