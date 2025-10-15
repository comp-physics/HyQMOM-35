# Examples - 3D Crossing Jets Simulations

All examples support **full parameter control** via both code defaults and command-line arguments.

## Available Examples

### 🎬 `run_3d_jets_timeseries.jl` - Interactive Time-Series (Recommended)
**Main example with interactive 3D visualization over time. Works in both serial and MPI parallel modes.**

```bash
# Serial run
julia --project=. examples/run_3d_jets_timeseries.jl

# Serial with custom parameters
julia --project=. examples/run_3d_jets_timeseries.jl --Np 60 --tmax 0.1 --snapshot-interval 5

# MPI parallel (4 ranks)
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100

# MPI parallel with custom parameters
mpiexec -n 8 julia --project=. examples/run_3d_jets_timeseries.jl --Np 120 --Nz 60 --snapshot-interval 5
```

**Features:**
- Automatic serial/MPI support (no separate MPI example needed!)
- Interactive GLMakie 3D viewer
- Time-series animation with play/pause
- Multiple quantities (density, U/V/W velocities)
- Isosurface visualization with positive/negative values
- Full parameter control via command-line

### 📊 `run_3d_crossing_jets.jl` - Static Plots
**Static matplotlib/PyPlot visualization (no GLMakie required). Also supports serial and MPI.**

```bash
# Serial
julia --project=. examples/run_3d_crossing_jets.jl
julia --project=. examples/run_3d_crossing_jets.jl --Np 60 --tmax 0.05

# MPI parallel
mpiexec -n 4 julia --project=. examples/run_3d_crossing_jets.jl --Np 100
```

**Features:**
- Automatic serial/MPI support
- Static 2D/3D plots (PyPlot/matplotlib)
- Multiple slice visualizations
- Centerline profiles
- Good for publications/reports

## Quick Start

### First Time Setup
```bash
cd HyQMOM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run with Defaults
```bash
julia --project=. examples/run_3d_jets_timeseries.jl
```

### Get Help
```bash
julia --project=. examples/run_3d_jets_timeseries.jl --help
```

## Parameter Control

### Command-Line Override System
**Priority (highest to lowest):**
1. Command-line arguments (`--Np 60`)
2. Code defaults in the example file
3. Global defaults from `parse_params.jl`

### Common Parameters

#### Grid & Domain
```bash
--Np N          # Grid resolution in x,y (default: 40)
--Nz N          # Grid resolution in z (default: 20)
--xmin X        # Domain x minimum (default: -0.5)
--xmax X        # Domain x maximum (default: 0.5)
# Similar for --ymin, --ymax, --zmin, --zmax
```

**Note:** Grid spacing `dx`, `dy`, `dz` are computed automatically!

#### Time & Physics
```bash
--tmax T        # Maximum time (default: 0.05)
--CFL C         # CFL number (default: 0.7)
--Ma M          # Mach number (default: 1.0)
--Kn K          # Knudsen number (default: 1.0)
```

#### Snapshots & Diagnostics
```bash
--snapshot-interval N    # Save every N steps (0=disabled, default varies)
--homogeneous-z BOOL     # Z-homogeneous mode (default: false)
```

### Usage Examples

#### Quick Low-Resolution Test
```bash
julia run_3d_jets_timeseries.jl --Np 20 --Nz 10 --tmax 0.01 --snapshot-interval 5
```

#### High-Resolution Run
```bash
mpiexec -n 4 julia run_3d_jets_mpi.jl --Np 80 --Nz 40 --tmax 0.2
```

#### Custom Domain
```bash
julia run_3d_jets_timeseries.jl \
  --xmin 0 --xmax 2 --ymin 0 --ymax 2 --zmin 0 --zmax 2 --Np 80
```

#### Parameter Sweep
```bash
for Ma in 0.5 0.7 1.0 1.5; do
  julia run_3d_jets_timeseries.jl --Ma $Ma --tmax 0.1 --snapshot-interval 5
done
```

#### Production Run (No Visualization)
```bash
julia run_3d_jets_timeseries.jl --Np 100 --tmax 1.0 --snapshot-interval 0
```

## Implementation Details

### Parameter Parsing System
Located in `parse_params.jl`:

```julia
include("parse_params.jl")

# Parse with code defaults + command-line overrides
params = parse_simulation_params(
    Np = 40,    # Code default, override with --Np N
    Nz = 20,    # Code default, override with --Nz N
    tmax = 0.05 # Code default, override with --tmax T
)

# Print parameter summary
print_params_summary(params, rank=rank)
```

### Key Functions
- `parse_simulation_params()`: Merges defaults + code + CLI
- `print_params_summary()`: Pretty-prints all parameters
- `get_default_params()`: Returns global defaults

## Tips & Best Practices

### For Quick Tests
- Use low resolution: `--Np 20 --Nz 10`
- Short time: `--tmax 0.01`
- Fewer snapshots: `--snapshot-interval 10`

### For Production
- Use MPI: `mpiexec -n 8 julia run_3d_jets_mpi.jl`
- Higher resolution: `--Np 100 --Nz 50`
- Disable snapshots if not needed: `--snapshot-interval 0`

### For Visualization
- Balance snapshot interval with memory
- Too many snapshots (interval=1) → uses lots of RAM
- Too few snapshots (interval=20) → choppy animation
- Sweet spot: `--snapshot-interval 2` to `5`

### For Parameter Sweeps
- Use command-line args in loops
- Save results to different directories
- Consider using `--snapshot-interval 0` for sweep

## Troubleshooting

### Out of Memory
```bash
# Reduce resolution
julia ... --Np 30 --Nz 15

# Increase snapshot interval
julia ... --snapshot-interval 10

# Disable snapshots
julia ... --snapshot-interval 0
```

### Simulation Crashes (NaN)
```bash
# Reduce Mach number
julia ... --Ma 0.7

# Reduce CFL
julia ... --CFL 0.5

# Increase resolution
julia ... --Np 60 --Nz 30
```

### GLMakie Viewer Doesn't Open
- Check GLMakie is installed: `using Pkg; Pkg.add("GLMakie")`
- Use static plots instead: `run_3d_crossing_jets.jl`
- Run on rank 0 only (already handled in MPI examples)

### Performance Issues
- Use MPI for large grids: `mpiexec -n 4 julia run_3d_jets_timeseries.jl --Np 100`
- Reduce snapshot frequency: `--snapshot-interval 10`
- Disable diagnostics: `--debug-output false`

## File Structure

```
examples/
├── README.md                    # This file
├── parse_params.jl              # Parameter parsing utilities (MPI-aware)
├── run_3d_jets_timeseries.jl   # Main: interactive time-series (serial/MPI)
└── run_3d_crossing_jets.jl     # Static PyPlot visualization (serial/MPI)
```

## Advanced Usage

### Creating Your Own Example

```julia
using HyQMOM
using MPI

include("parse_params.jl")

MPI.Init()

# Parse parameters with your defaults
params = parse_simulation_params(
    Np = 40,
    # ... your custom defaults
)

print_params_summary(params, rank=MPI.Comm_rank(MPI.COMM_WORLD))

# Run simulation
if params.snapshot_interval > 0
    snapshots, grid = run_simulation_with_snapshots(params; 
                                                     snapshot_interval=params.snapshot_interval)
else
    M_final, final_time, time_steps, grid = simulation_runner(params)
end

MPI.Finalize()
```

### Modifying Defaults

Edit `parse_params.jl` function `get_default_params()` to change global defaults.

### Adding New Parameters

1. Add to `get_default_params()` in `parse_params.jl`
2. Add to help text in `print_help()`
3. Use in your example with `params.your_new_param`

## References

- Main documentation: `../README.md`
- Parameter details: `parse_params.jl` source
- Grid spacing removal: See commit history (dx/dy/dz now automatic)

## Questions?

Run any example with `--help` to see all available options:
```bash
julia --project=. examples/run_3d_jets_timeseries.jl --help
```
