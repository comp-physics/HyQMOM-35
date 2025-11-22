# HyQMOM.jl User Guide

A complete guide for running the HyQMOM Julia implementation, designed for users new to Julia.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Installing Julia](#installing-julia)
3. [Getting the Code](#getting-the-code)
4. [Setting Up the Environment](#setting-up-the-environment)
5. [Running Your First Simulation](#running-your-first-simulation)
6. [Understanding the Output](#understanding-the-output)
7. [Customizing Parameters](#customizing-parameters)
8. [Running with MPI (Parallel)](#running-with-mpi-parallel)
9. [Visualization Options](#visualization-options)
10. [Troubleshooting](#troubleshooting)
11. [Advanced Usage](#advanced-usage)

---

## Prerequisites

### What You Need

- **Operating System**: Linux, macOS, or Windows
- **Memory**: At least 4 GB RAM (8 GB+ recommended for larger simulations)
- **Disk Space**: About 2 GB for Julia and packages
- **Time**: 15-30 minutes for first-time setup

### What is Julia?

Julia is a modern programming language designed for scientific computing. It's:
- Fast (comparable to C/Fortran)
- Easy to use (like Python or MATLAB)
- Free and open-source

You don't need to know Julia to run these simulations! This guide walks you through everything.

---

## Installing Julia

### Step 1: Download Julia

1. Go to [https://julialang.org/downloads/](https://julialang.org/downloads/)
2. Download Julia **1.9 or later** (1.10 is recommended)
3. Choose the installer for your operating system:
   - **Windows**: Download the 64-bit `.exe` installer
   - **macOS**: Download the `.dmg` file
   - **Linux**: Download the `.tar.gz` file

### Step 2: Install Julia

**On Windows:**
1. Run the downloaded `.exe` file
2. Follow the installation wizard
3. Check "Add Julia to PATH" during installation
4. Click "Install"

**On macOS:**
1. Open the downloaded `.dmg` file
2. Drag Julia to your Applications folder
3. Open Terminal and add Julia to your PATH:
   ```bash
   sudo ln -s /Applications/Julia-1.10.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
   ```
   (Adjust version number if different)

**On Linux:**
1. Extract the downloaded file:
   ```bash
   tar -xvzf julia-1.10.0-linux-x86_64.tar.gz
   ```
2. Move to `/opt/`:
   ```bash
   sudo mv julia-1.10.0 /opt/
   ```
3. Create symbolic link:
   ```bash
   sudo ln -s /opt/julia-1.10.0/bin/julia /usr/local/bin/julia
   ```

### Step 3: Verify Installation

Open a terminal (or Command Prompt on Windows) and type:

```bash
julia --version
```

You should see something like:
```
julia version 1.10.0
```

If you see this, Julia is installed correctly!

---

## Getting the Code

### Step 1: Clone the Repository

If you have Git installed:

```bash
cd ~  # or wherever you want to put the code
git clone https://github.com/comp-physics/HyQMOM.jl.git
cd HyQMOM.jl
```

If you don't have Git:
1. Go to the GitHub repository page
2. Click the green "Code" button
3. Click "Download ZIP"
4. Extract the ZIP file
5. Open a terminal in the extracted folder

### Step 2: You're Ready!

The repository root is the Julia package directory, containing all the source code, tests, examples, and documentation.

---

## Setting Up the Environment

### Step 1: Install Package Dependencies

This is a **one-time setup**.
In the terminal, from the repository root directory, run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

**What this does:**
- `julia`: Starts Julia
- `--project=.`: Uses the current directory's project settings
- `-e`: Executes the command that follows
- `Pkg.instantiate()`: Installs all required packages

**This will take 5-15 minutes** the first time, as it downloads and compiles packages.

You'll see output like:
```
Updating registry at `~/.julia/registries/General.toml`
Installed MPI ──────────── v0.20.23
Installed GLMakie ──────── v0.10.18
...
Precompiling project...
```

When it finishes, you'll see:
```
✓ All packages installed successfully
```

### Step 2: Verify Setup

Test that everything works:

```bash
julia --project=. -e 'using HyQMOM; println("Setup complete!")'
```

If you see `Setup complete!`, you're ready to run simulations!

---

## Running Your First Simulation

### Quick Test Run

Let's run a small, fast simulation to make sure everything works.

**Command:**

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 20 --Ny 20 --tmax 0.01
```

**What this command means:**

- `julia --project=.`: Run Julia with this project's settings
- `examples/run_3d_jets_timeseries.jl`: The simulation script to run
- `--Nx 20 --Ny 20`: Use a 20×20 grid in x-y plane (small for speed)
- `--tmax 0.01`: Run until time t=0.01 (very short)

**What you'll see:**

```
Starting time evolution...
  Grid: 20x20x20, Ranks: 1, Local: 20x20x20
  tmax: 0.01, CFL: 0.7, Ma: 1.0, Kn: 1.0
Step    1: t = 0.0050, dt = 5.0000e-03, wall = 0.2 s
Step    2: t = 0.0100, dt = 5.0000e-03, wall = 0.2 s
Time evolution complete: 2 steps, t = 0.01
```

**If you see a 3D visualization window open:** Congratulations! The simulation ran successfully. You can:
- Drag with mouse to rotate
- Scroll to zoom
- Use the slider to see time evolution
- Click buttons to view different quantities (density, velocity, etc.)

**If no window opens:** That's okay! The simulation still ran. See [Visualization Options](#visualization-options) below.

### Understanding the Simulation

This simulation models **two jets of gas crossing in 3D space**:
- The jets collide and interact
- The code tracks 35 moments of the 3D velocity distribution
- It solves the flow equations with BGK collision operator and closure via HyQMOM

---

## Understanding the Output

### Terminal Output Explained

```
Starting time evolution...
  Grid: 20x20x20, Ranks: 1, Local: 20x20x20
```
- **Grid**: The spatial resolution (20 points in each direction)
- **Ranks**: Number of parallel processes (1 = serial)
- **Local**: Grid size per process

```
  tmax: 0.01, CFL: 0.7, Ma: 1.0, Kn: 1.0
```
- **tmax**: Maximum simulation time
- **CFL**: Stability parameter (smaller = more stable)
- **Ma**: Mach number (flow speed / speed of sound)
- **Kn**: Knudsen number (mean free path / length scale)

```
Step    1: t = 0.0050, dt = 5.0000e-03, wall = 0.2 s
```
- **Step**: Time step number
- **t**: Current simulation time
- **dt**: Time step size (automatically adjusted)
- **wall**: Real time elapsed

### Files Created

After running, you may see:
- `results.jld2`: Saved simulation data (if snapshots enabled)
- Console output showing convergence and diagnostics

---

## Customizing Parameters

### Basic Parameters

You can customize the simulation with command-line options. Here are the most important ones:

#### Grid Resolution

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40 --Nz 40
```

- `--Nx 40`: Grid size in x direction (default: 40)
- `--Ny 40`: Grid size in y direction (default: 40)
- `--Nz 40`: Grid size in z direction (default: 20)
- **Larger = more accurate but slower**
- Memory use: ~4 GB for 60×60×40, ~16 GB for 100×100×60

#### Simulation Time

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --tmax 0.2
```

- `--tmax 0.2`: Run until t=0.2 (default: 0.05)
- **Longer = more evolution but takes more time**

#### Physical Parameters

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Ma 0.5 --Kn 2.0
```

- `--Ma 0.5`: Mach number (jet speed / sound speed)
  - Lower = subsonic, Higher = supersonic
- `--Kn 2.0`: Knudsen number (rarefaction)
  - Lower = continuum, Higher = free molecular flow

#### Stability

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --CFL 0.5
```

- `--CFL 0.5`: CFL stability parameter (default: 0.7)
- **If simulation crashes, try lower CFL (e.g., 0.3)**

### Example Combinations

**Quick test (2-3 minutes):**
```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 20 --Ny 20 --tmax 0.01 --snapshot-interval 1
```

**Standard run (15-30 minutes):**
```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40 --tmax 0.1
```

**High quality (1-2 hours):**
```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60 --tmax 0.2 --CFL 0.6
```

**Production (several hours):**
```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 100 --Ny 100 --tmax 0.5 --CFL 0.5
```

### Getting Help

To see all available options:

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --help
```

---

## Running with MPI (Parallel)

MPI allows the simulation to run on multiple CPU cores simultaneously, making it much faster.

### Step 1: Install MPI

**On macOS (with Homebrew):**
```bash
brew install mpich
```

**On Ubuntu/Debian Linux:**
```bash
sudo apt-get install mpich libmpich-dev
```

**On Windows:**
- Download Microsoft MPI from [Microsoft's website](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi)
- Install both `msmpisetup.exe` and `msmpisdk.msi`

### Step 2: Configure MPI.jl

```bash
julia --project=. -e 'using Pkg; Pkg.build("MPI")'
```

### Step 3: Run in Parallel

Use `mpiexec` to run with multiple processes:

```bash
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60
```

**What this means:**
- `mpiexec -n 4`: Run with 4 parallel processes
- The grid is divided among processes automatically
- **Each process needs memory**, so don't use too many

**How many processes to use?**
- **Start with 2-4** to test
- Good rule: Number of CPU cores / 2
- Example: 8-core CPU → try `-n 4`
- Don't exceed your number of CPU cores

**Expected speedup:**
- 2 processes: ~1.7× faster
- 4 processes: ~3× faster
- 8 processes: ~5× faster

### Parallel Examples

**Medium simulation, 4 cores (10-15 minutes):**
```bash
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60 --tmax 0.1
```

**Large simulation, 8 cores (30-60 minutes):**
```bash
mpiexec -n 8 julia --project=. examples/run_3d_jets_timeseries.jl --Nx 80 --Ny 80 --tmax 0.2
```

---

## Visualization Options

### Interactive 3D Visualization

**Default behavior:** The simulation opens an interactive 3D window showing:
- Isosurfaces of density, velocity, pressure, temperature
- Time slider to see evolution
- Play/pause controls
- Rotation and zoom with mouse

**Controls:**
- **Left mouse drag**: Rotate view
- **Scroll wheel**: Zoom in/out
- **Time slider**: Move through simulation snapshots
- **Play button**: Animate time evolution
- **Quantity buttons**: Switch between density, U, V, W velocities, etc.

### Disabling Visualization

For faster runs or headless systems (no display):

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40 --snapshot-interval 0
```

- `--snapshot-interval 0`: Disables snapshot collection (no visualization)
- Simulation runs faster without storing intermediate states

Or set environment variable:
```bash
export HYQMOM_SKIP_PLOTTING=true
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40
```

### Saving Snapshots

To save simulation data for later analysis:

```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40 --snapshot-interval 5
```

- `--snapshot-interval 5`: Save every 5 time steps
- Creates `results.jld2` file with simulation data
- Can load later for analysis or visualization


## Troubleshooting

### Problem: "Package not found" or "Module not found"

**Solution:**
```bash
cd HyQMOM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Make sure you're in the `HyQMOM.jl` directory and using `--project=.`

### Problem: "No display found" or GLMakie errors

This happens when running on a system without a graphics display (like a server).

**Solution:** Disable visualization:
```bash
export HYQMOM_SKIP_PLOTTING=true
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40
```

Or use snapshot mode:
```bash
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40 --snapshot-interval 0
```

### Problem: Simulation crashes with NaN values

This means numerical instability.

**Solutions:**
1. **Reduce CFL number:**
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --CFL 0.3
   ```

2. **Increase grid resolution:**
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60
   ```

3. **Reduce Mach number:**
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --Ma 0.5
   ```

### Problem: "Out of memory" errors

The simulation needs too much RAM.

**Solutions:**
1. **Reduce grid size:**
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --Nx 30 --Ny 30
   ```

2. **Disable snapshots:**
   ```bash
   julia --project=. examples/run_3d_jets_timeseries.jl --snapshot-interval 0
   ```

3. **Use MPI to distribute memory:**
   ```bash
   mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40
   ```

### Problem: MPI doesn't work

**Check MPI installation:**
```bash
mpiexec --version
```

If not found, install MPI (see [Running with MPI](#running-with-mpi-parallel)).

**Rebuild MPI.jl:**
```bash
julia --project=. -e 'using Pkg; Pkg.build("MPI")'
```

### Problem: Simulation is too slow

**Speed it up:**
1. **Use smaller grid:** `--Nx 30 --Ny 30` instead of `--Nx 60 --Ny 60`
2. **Shorter time:** `--tmax 0.05` instead of `--tmax 0.2`
3. **Use MPI:** `mpiexec -n 4 julia ...`
4. **Disable snapshots:** `--snapshot-interval 0`

### Problem: "Command not found: julia"

Julia is not in your system PATH.

**Find Julia:**
```bash
# On macOS/Linux
which julia

# On Windows
where julia
```

If not found, add Julia to PATH:
- **macOS/Linux:** Edit `~/.bashrc` or `~/.zshrc` and add:
  ```bash
  export PATH="/path/to/julia/bin:$PATH"
  ```
- **Windows:** Add Julia's `bin` directory to system PATH in Environment Variables

---

## Advanced Usage

### Custom Initial Conditions

You can create custom jet configurations:

```bash
julia --project=. examples/run_3d_custom_jets.jl --config triple-jet
```

Available configurations:
- `crossing`: Two perpendicular jets (default)
- `crossing`: Two perpendicular jets (default)
- `triple-jet`: Three jets at 120° angles
- `quad-jet`: Four jets converging to center
- `vertical-jet`: Single vertical jet
- `spiral`: Spiral jet pattern

### Running Tests

To verify the installation:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This runs the test suite (~2 minutes).

### Checking Against MATLAB Results

To validate the Julia code against MATLAB reference:

```bash
julia --project=. test/test_golden_files.jl
```

This compares Julia results to MATLAB golden files.

### Using from Julia Scripts

You can write your own Julia scripts:

```julia
using HyQMOM

# Set parameters
params = (
    Nx = 40,
    Ny = 40,
    Nz = 40,
    tmax = 0.1,
    Ma = 1.0,
    Kn = 1.0,
    CFL = 0.7,
)

# Run simulation
M_final, final_time, steps, grid = simulation_runner(params)

println("Simulation complete!")
println("Final time: $final_time")
println("Time steps: $steps")
```

Save as `my_simulation.jl` and run:
```bash
julia --project=. my_simulation.jl
```

### Analyzing Results

If you saved snapshots, load them in Julia:

```julia
using JLD2

# Load saved data
data = load("results.jld2")
snapshots = data["snapshots"]
grid = data["grid"]

# Access specific snapshot
snapshot_10 = snapshots[10]
density = snapshot_10[:, :, :, 1]  # First moment is density

println("Max density: ", maximum(density))
```

---

## Summary: Quick Reference

### First Time Setup (Once)

```bash
# 1. Install Julia from julialang.org
# 2. Clone repository
git clone https://github.com/comp-physics/HyQMOM.jl.git
cd HyQMOM.jl

# 3. Install packages (takes 10-15 min)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Every Time You Run

```bash
# Navigate to repository root
cd path/to/HyQMOM.jl

# Run simulation
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40
```

### Common Commands

```bash
# Quick test (2 min)
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 20 --Ny 20 --tmax 0.01

# Standard run (15 min)
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40 --tmax 0.1

# High quality (1 hour)
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60 --tmax 0.2

# Parallel (4 cores)
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60

# Without visualization
julia --project=. examples/run_3d_jets_timeseries.jl --Nx 40 --Ny 40 --snapshot-interval 0
```

### Getting Help

```bash
# See all options
julia --project=. examples/run_3d_jets_timeseries.jl --help

# Check installation
julia --project=. -e 'using HyQMOM; println("OK")'

# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

---

## Further Resources

- **Julia Documentation**: [https://docs.julialang.org](https://docs.julialang.org)
- **Main README**: See `README.md` for technical details
- **Example Documentation**: See `HyQMOM.jl/examples/README.md` for more examples

---

## Tips for Success

1. **Start small**: Use `--Nx 20 --Ny 20 --tmax 0.01` for your first runs
2. **Check memory**: Monitor RAM usage with `htop` (Linux/Mac) or Task Manager (Windows)
3. **Save your work**: Use descriptive output names and save parameter combinations
4. **Be patient**: First run is slower due to compilation (subsequent runs are faster)
5. **Use MPI**: Parallel runs are much faster for large problems
6. **Read errors**: Error messages usually tell you exactly what's wrong
7. **Keep backups**: Save `results.jld2` files from important runs
