# Interactive 3D Visualization Tutorial

This tutorial walks you through HyQMOM.jl's interactive 3D visualization capabilities, from basic usage to advanced customization.

## Overview

HyQMOM.jl provides real-time 3D visualization using GLMakie, allowing you to:
- Visualize 3D isosurfaces of physical quantities (density, velocity, pressure, temperature)
- Animate time evolution with interactive controls
- Explore simulation results with mouse-driven camera controls
- Switch between different physical quantities in real-time

## Basic Usage

### Running Your First Visualization

Start with the interactive time-series example:

```bash
julia --project=. examples/run_3d_jets_timeseries.jl
```

This will:
1. Run a 3D crossing jets simulation
2. Automatically launch the GLMakie visualization window
3. Display density isosurfaces by default
4. Save snapshots to a `.jld2` file

### Understanding the Interface

The visualization window contains several interactive elements:

**Time Controls:**
- **Time slider**: Drag to navigate through simulation snapshots
- **Play button**: Start/stop time animation
- **Reset button**: Return to the first timestep

**Quantity Selection:**
- **Density**: Shows mass density isosurfaces
- **U velocity**: X-component of velocity
- **V velocity**: Y-component of velocity  
- **W velocity**: Z-component of velocity
- **Pressure**: Pressure field
- **Temperature**: Temperature field

**Visualization Controls:**
- **Isosurface sliders**: Adjust the threshold values for isosurface rendering
- **Mouse controls**: Left-drag to rotate, scroll to zoom, right-drag to pan

## Step-by-Step Walkthrough

### Step 1: Launch a Quick Simulation

For this tutorial, we'll use a fast, low-resolution simulation:

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Np 30 --tmax 0.05 --snapshot-interval 2
```

**Parameters explained:**
- `--Np 30`: Use 30×30 grid points (faster computation)
- `--tmax 0.05`: Run for shorter time (faster completion)
- `--snapshot-interval 2`: Save every 2nd timestep (more animation frames)

### Step 2: Explore the Density Field

When the visualization opens:

1. **Observe the initial state**: You should see two jets crossing in 3D space
2. **Press Play**: Watch the jets evolve and interact over time
3. **Adjust the time slider**: Manually step through different time points
4. **Modify isosurface levels**: Use the sliders to show more or less of the density field

### Step 3: Switch to Velocity Visualization

1. **Click "U velocity"**: Switch to x-component of velocity
2. **Observe the differences**: Velocity fields often show different structures than density
3. **Try other components**: Click "V velocity" and "W velocity" to see the full 3D flow
4. **Compare with density**: Switch back and forth to understand the relationship

### Step 4: Camera Navigation

Master the 3D camera controls:

1. **Rotation**: Left-click and drag to rotate around the domain
2. **Zooming**: Scroll wheel to zoom in/out
3. **Panning**: Right-click and drag to pan the view
4. **Find good angles**: Rotate to see the jet interaction from different perspectives

## Advanced Features

### Customizing Isosurface Levels

The isosurface sliders control which parts of the field are visible:

- **Lower threshold**: Shows more of the field (including low-intensity regions)
- **Higher threshold**: Shows only high-intensity regions (cleaner visualization)
- **Multiple levels**: Some quantities support multiple isosurface levels simultaneously

**Tips:**
- For density: Try values between 0.5 and 2.0 times the background density
- For velocity: Experiment with 10-50% of the maximum velocity magnitude
- For pressure/temperature: Use values around the mean ± 1-2 standard deviations

### Working with Saved Data

After running a simulation, you can re-visualize the results:

```bash
# The simulation saves a .jld2 file automatically
julia visualize_jld2.jl

# Or specify the file directly
julia visualize_jld2.jl snapshots_crossing_Ma1.0_t0.05_N30.jld2
```

This is useful for:
- Re-examining results without re-running simulations
- Sharing visualizations with colleagues
- Creating animations or screenshots for presentations

### Programmatic Visualization

For advanced users, you can control visualization programmatically:

```julia
using HyQMOM, JLD2, GLMakie

# Load simulation data
@load "snapshots_file.jld2" snapshots grid params params_with_ic

# Launch interactive viewer
interactive_3d_timeseries_streaming("snapshots_file.jld2", grid, params_with_ic)

# Or visualize a specific timestep
final_state = snapshots[end]
interactive_standardized_scatter(final_state, grid)
```

## Parallel Visualization

### MPI Parallel Simulations

When running with MPI, visualization automatically appears on rank 0:

```bash
# Run parallel simulation - visualization will appear on rank 0
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 60
```

**Key points:**
- Only rank 0 opens the visualization window
- All ranks participate in the computation
- Snapshots are saved from the combined parallel data
- No special configuration needed

### Headless Systems

For systems without display capabilities (HPC clusters, remote servers):

```bash
# Skip visualization during simulation
julia --project=. examples/run_3d_jets_timeseries.jl --no-viz true

# Or use environment variable
export HYQMOM_SKIP_PLOTTING=true
julia --project=. examples/run_3d_jets_timeseries.jl

# Later, on a system with display:
julia visualize_jld2.jl snapshots_file.jld2
```

## Troubleshooting

### Common Issues

**GLMakie window doesn't appear:**
```bash
# Test GLMakie installation
julia --project=. -e 'using GLMakie; GLMakie.activate!()'

# On remote systems, enable X11 forwarding
ssh -Y username@hostname

# Check for display issues
echo $DISPLAY
```

**Slow visualization performance:**
```bash
# Reduce resolution
julia ... --Np 20

# Increase snapshot interval (fewer frames)
julia ... --snapshot-interval 5

# Use software rendering if GPU issues
export LIBGL_ALWAYS_SOFTWARE=1
```

**Visualization crashes or freezes:**
```bash
# Try software rendering
export LIBGL_ALWAYS_SOFTWARE=1

# Reduce memory usage
julia ... --Np 30 --snapshot-interval 10

# Check system resources
htop  # Monitor CPU and memory usage
```

### Platform-Specific Notes

**macOS:**
- GLMakie should work out of the box
- If issues occur, try updating graphics drivers

**Linux:**
- Ensure OpenGL drivers are installed
- For headless systems, use `LIBGL_ALWAYS_SOFTWARE=1`
- X11 forwarding: `ssh -Y` or `ssh -X`

**Windows:**
- GLMakie typically works with modern graphics drivers
- WSL users may need X11 server (VcXsrv, Xming)

## Tips and Best Practices

### Effective Visualization Workflow

1. **Start small**: Use low resolution (`--Np 20-30`) for initial exploration
2. **Iterate parameters**: Adjust physics parameters and observe changes
3. **Save interesting cases**: Keep `.jld2` files for later analysis
4. **Document findings**: Take screenshots or notes of interesting phenomena

### Creating Animations

For creating movies or GIFs:

1. **Run simulation with frequent snapshots**: `--snapshot-interval 1`
2. **Use consistent parameters**: Keep resolution and time steps regular
3. **Post-process**: Use the saved `.jld2` files to create custom animations

### Performance Optimization

- **Resolution vs. quality**: Balance computational cost with visualization detail
- **Snapshot frequency**: More snapshots = smoother animation but larger files
- **Isosurface complexity**: Simpler isosurfaces render faster

## Next Steps

After mastering basic visualization:

1. **Explore different initial conditions**: Try the custom jets examples
2. **Experiment with physics parameters**: Vary Mach and Knudsen numbers
3. **Scale up**: Run higher resolution simulations with MPI
4. **Customize**: Modify visualization code for specific analysis needs

See the User Guide for more simulation options and the API Reference for programmatic control.
