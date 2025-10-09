# HyQMOM.jl Examples

This directory contains example scripts demonstrating how to use the HyQMOM.jl package.

## Quick Start

### Simple Example
Run a small 40x40 grid simulation:

```bash
cd examples
julia --project=.. run_simple.jl
```

With MPI (2 ranks):
```bash
mpiexec -n 2 julia --project=.. run_simple.jl
```

### Visualization Example
Run a simulation with PyPlot visualization (equivalent to MATLAB's plotting):

```bash
cd examples
julia --project=.. run_with_plots.jl
```

This example generates all final result plots (Figures 2-12):
- Figure 2: Moment line plots along diagonal
- Figure 3: Central moment line plots
- Figure 4: Standardized moment line plots
- Figure 9: Contour plots (12 panels)
- Figure 10: C-moment contour plots (16 panels)
- Figure 11: S-moment contour plots (12 panels)
- Figure 12: Hyperbolicity plots (9 panels)

**Requirements**: PyPlot must be installed (`using Pkg; Pkg.add("PyPlot")`)

## Using the Top-Level main.jl

The repository includes a convenient top-level `main.jl` file (similar to MATLAB's `main.m`) that you can run from the repository root:

### Single Rank
```bash
# From repository root
julia main.jl

# With custom parameters
julia main.jl --Np 40 --tmax 0.1
```

### MPI Parallel
```bash
# 2 ranks
mpiexec -n 2 julia main.jl

# 4 ranks with custom parameters
mpiexec -n 4 julia main.jl --Np 80 --tmax 0.1 --Ma 2.0

# 8 ranks for production run
mpiexec -n 8 julia main.jl --Np 240 --tmax 0.04 --Kn 0.01
```

### Match MATLAB Golden File Parameters
To run with the same parameters as the MATLAB golden file (for comparison):

```bash
julia main.jl --Np 20 --tmax 0.1 --Kn 1.0 --Ma 0.0 --CFL 0.5
```

## Command-Line Arguments

The `main.jl` script accepts the following arguments:

- `--Np <int>`: Grid size (default: 120)
- `--tmax <float>`: Maximum simulation time (default: 0.02)
- `--Ma <float>`: Mach number (default: 2.0)
- `--Kn <float>`: Knudsen number (default: 0.01)
- `--CFL <float>`: CFL number (default: 0.9)
- `--flag2D <int>`: 2D flag (default: 0)
- `--output <string>`: Output file name (default: "results.jld2")
- `--help`: Show help message

## Output

Simulation results are saved in JLD2 format containing:
- `M`: Final moment field (Np x Np x 35)
- `t`: Final time
- `steps`: Number of time steps taken
- `grid`: Grid information

## Tips

1. **Start small**: Test with `--Np 40 --tmax 0.1` before running large simulations
2. **MPI scaling**: For best performance, use powers of 2 for the number of ranks (2, 4, 8, etc.)
3. **Grid decomposition**: Each rank needs at least 10x10 interior points
4. **Memory**: Large grids (Np > 200) may require significant RAM

## Comparison with MATLAB

| MATLAB | Julia |
|--------|-------|
| `main.m` | `main.jl` (repository root) |
| `matlab main.m 40 0.1` | `julia main.jl --Np 40 --tmax 0.1` |
| `parpool('local', 4)` | `mpiexec -n 4 julia main.jl` |

