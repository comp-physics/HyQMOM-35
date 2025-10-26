# User Guide

This comprehensive guide covers all aspects of using HyQMOM.jl for computational fluid dynamics simulations.

## Command-Line Parameters

HyQMOM.jl examples support extensive command-line configuration. All parameters have sensible defaults but can be customized for your specific needs.

### Grid and Time Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--Nx N` | Grid resolution in x direction | 40 | `--Nx 100` |
| `--Ny N` | Grid resolution in y direction | 40 | `--Ny 100` |
| `--Nz N` | Grid resolution in z direction | 20 | `--Nz 60` |
| `--tmax T` | Maximum simulation time | 0.05 | `--tmax 0.5` |
| `--CFL C` | CFL number for stability | 0.7 | `--CFL 0.5` |

### Physics Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--Ma M` | Mach number (jet velocity / thermal velocity) | 0.0 | `--Ma 1.5` |
| `--Kn K` | Knudsen number (mean free path / characteristic length) | 1.0 | `--Kn 0.5` |

**Physical Regimes:**
- **Low Kn (< 0.1)**: Continuum/hydrodynamic regime with frequent collisions
- **Moderate Kn (0.1-10)**: Transitional regime with mixed kinetic-continuum behavior  
- **High Kn (> 10)**: Free molecular/kinetic regime with rare collisions
- **Low Ma (< 1)**: Subsonic flows, moments remain close to equilibrium
- **High Ma (> 1)**: Supersonic flows, strong non-equilibrium effects, crossing jets

### Visualization and Output

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--snapshot-interval N` | Save every N steps (0=disabled) | 0 | `--snapshot-interval 5` |
| `--no-viz true/false` | Skip interactive visualization | false | `--no-viz true` |

### Getting Help

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --help
```

## Simulation Workflows

### Interactive Development Workflow

1. **Quick test**: Start with low resolution for rapid iteration
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --Nx 20 --Ny 20 --tmax 0.01
   ```

2. **Parameter exploration**: Adjust physics parameters
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --Ma 0.5 --Kn 2.0
   ```

3. **Production run**: Scale up resolution and time
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --Nx 80 --Ny 80 --tmax 0.5
   ```

### Batch/HPC Workflow

1. **Prepare parameters**: Use `--no-viz true` for headless systems
2. **Submit job**: Use MPI for parallel execution
3. **Post-process**: Transfer `.jld2` files for visualization

```bash
# On HPC system
export HYQMOM_SKIP_PLOTTING=true
mpiexec -n 16 julia --project=. examples/run_3d_jets_timeseries.jl \
  --Nx 200 --Ny 200 --tmax 1.0 --no-viz true --snapshot-interval 5

# On workstation
scp cluster:path/to/snapshots_*.jld2 .
julia visualize_jld2.jl
```

## Visualization Modes

### Interactive 3D Visualization

The primary visualization mode uses GLMakie for real-time 3D rendering:

**Features:**
- Real-time 3D isosurface visualization
- Time slider with play/pause animation
- Multiple quantities (density, U/V/W velocities, pressure, temperature)
- Interactive camera controls (rotate, zoom, pan)
- Adjustable isosurface levels

**Controls:**
- **Time slider**: Navigate through simulation snapshots
- **Play button**: Animate time evolution
- **Quantity buttons**: Switch between physical fields
- **Isosurface sliders**: Adjust visualization thresholds
- **Mouse**: Left-drag to rotate, scroll to zoom, right-drag to pan

### Static Plots

For publication-quality figures, use the static plotting examples:

```bash
julia --project=. examples/run_3d_crossing_jets.jl
```

This generates high-quality matplotlib figures suitable for papers and presentations.

### Post-Processing Visualization

The `visualize_jld2.jl` script provides flexible post-processing:

```bash
# Auto-detect files
julia visualize_jld2.jl

# Specify file
julia visualize_jld2.jl my_simulation.jld2

# Multiple files (will prompt for selection)
julia visualize_jld2.jl *.jld2
```

## Initial Conditions and Configurations

### Predefined Configurations

HyQMOM.jl includes several predefined initial condition configurations:

```bash
# Crossing jets (default)
julia --project=. examples/run_3d_custom_jets.jl --config crossing

# Triple jet configuration
julia --project=. examples/run_3d_custom_jets.jl --config triple-jet

# Four jets converging to center
julia --project=. examples/run_3d_custom_jets.jl --config quad-jet

# Vertical jet
julia --project=. examples/run_3d_custom_jets.jl --config vertical-jet

# Spiral configuration
julia --project=. examples/run_3d_custom_jets.jl --config spiral
```

### Custom Initial Conditions

You can define custom initial conditions by modifying the configuration system. See `examples/CUSTOM_INITIAL_CONDITIONS.md` for detailed instructions on:

- Defining new jet configurations
- Specifying center positions, widths, and velocities
- Creating complex multi-jet scenarios

## Performance Optimization

### Memory Management

**Reduce memory usage:**
```bash
# Lower resolution
--Nx 30 --Ny 30 --Nz 15

# Reduce snapshot frequency
--snapshot-interval 10

# Disable snapshots entirely
--snapshot-interval 0
```

**Monitor memory usage:**
- Use system monitoring tools (`htop`, `top`)
- Check MPI rank memory distribution
- Consider domain decomposition efficiency

### Numerical Stability

**For stable simulations:**
```bash
# Reduce CFL number
--CFL 0.5

# Lower Mach number
--Ma 0.7

# Increase resolution (better resolves gradients)
--Nx 60 --Ny 60 --Nz 30
```

**Signs of instability:**
- NaN values in output
- Exponentially growing fields
- Simulation crashes

### Parallel Efficiency

**Optimal MPI usage:**
- Use powers of 2 for number of ranks when possible
- Balance computation vs. communication overhead
- Monitor load balancing across ranks

```bash
# Good scaling examples
mpiexec -n 4 julia ... --Nx 80 --Ny 80    # 40x40 per rank (2x2 decomposition)
mpiexec -n 16 julia ... --Nx 120 --Ny 120  # 30x30 per rank (4x4 decomposition)
mpiexec -n 64 julia ... --Nx 160 --Ny 160  # 20x20 per rank (8x8 decomposition)
```

## Environment Variables

### Plotting Control

```bash
# Skip all plotting dependencies (useful for CI/HPC)
export HYQMOM_SKIP_PLOTTING=true

# Alternative: CI auto-detection
export CI=true
```

### Julia Performance

```bash
# Set number of threads
export JULIA_NUM_THREADS=4

# Optimize compilation
export JULIA_CPU_TARGET="native"
```

## File Formats and Data

### Snapshot Files (.jld2)

HyQMOM.jl saves simulation data in JLD2 format, which includes:

- **snapshots**: Array of solution states at different times
- **grid**: Spatial grid information
- **params**: Simulation parameters
- **params_with_ic**: Parameters including initial conditions

### Accessing Data Programmatically

```julia
using JLD2, HyQMOM

# Load simulation data
@load "snapshots_file.jld2" snapshots grid params params_with_ic

# Access specific timestep
final_state = snapshots[end]

# Extract physical quantities
density = final_state[:, :, :, 1]  # First moment (density)
u_velocity = final_state[:, :, :, 2] ./ density
v_velocity = final_state[:, :, :, 3] ./ density
w_velocity = final_state[:, :, :, 4] ./ density
```

## Troubleshooting

### Common Issues

**GLMakie visualization problems:**
- Ensure X11 forwarding: `ssh -Y user@host`
- Check OpenGL support: `julia -e 'using GLMakie; GLMakie.activate!()'`
- Try software rendering: `export LIBGL_ALWAYS_SOFTWARE=1`

**MPI execution issues:**
- Verify MPI installation: `mpiexec --version`
- Check MPI.jl configuration: `julia -e 'using MPI; println(MPI.versioninfo())'`
- Try fewer ranks: `mpiexec -n 2 julia ...`

**Memory problems:**
- Reduce resolution: `--Nx 30 --Ny 30`
- Increase snapshot interval: `--snapshot-interval 10`
- Use more MPI ranks to distribute memory

**Numerical instabilities:**
- Reduce CFL: `--CFL 0.5`
- Lower Mach number: `--Ma 0.7`
- Increase resolution: `--Nx 60 --Ny 60`

### Getting Help

1. Check the API Reference for function documentation
2. Review [example configurations](https://github.com/comp-physics/HyQMOM.jl/tree/main/examples)
3. Open an issue on [GitHub](https://github.com/comp-physics/HyQMOM.jl/issues)
