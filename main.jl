#!/usr/bin/env julia
"""
Main entry point for HyQMOM 3D simulation - Julia version

This is the Julia equivalent of main.m in MATLAB.
It provides a convenient top-level interface to run simulations.

# Usage

## Single rank (no MPI):
```bash
julia main.jl
julia main.jl --Np 40 --tmax 0.1
```

## MPI parallel:
```bash
mpiexec -n 2 julia main.jl
mpiexec -n 4 julia main.jl --Np 80 --tmax 0.1
mpiexec -n 8 julia main.jl --Np 240 --tmax 0.04 --Ma 2.0 --Kn 0.01
```

# Command-line Arguments
- `--Np`: Grid size (default: 120)
- `--tmax`: Maximum simulation time (default: 0.02)
- `--Ma`: Mach number (default: 2.0)
- `--Kn`: Knudsen number (default: 0.01)
- `--CFL`: CFL number (default: 0.9)
- `--flag2D`: 2D flag (default: 0)
- `--output`: Output file name (default: "results.jld2")
- `--help`: Show this help message

# Examples
```bash
# Default parameters (Np=120, tmax=0.02, Ma=2.0, Kn=0.01)
julia main.jl

# Small test case
julia main.jl --Np 40 --tmax 0.1

# Production run with MPI
mpiexec -n 8 julia main.jl --Np 240 --tmax 0.04

# Match MATLAB golden file parameters (1 rank, Np=20, tmax=0.1)
julia main.jl --Np 20 --tmax 0.1 --Kn 1.0 --Ma 0.0 --CFL 0.5
```

# Output
Results are saved to the specified output file (default: results.jld2)
containing:
- M: Final moment field
- t: Final time
- steps: Number of time steps
- grid: Grid information
"""

# Add the HyQMOM.jl package to the load path
push!(LOAD_PATH, joinpath(@__DIR__, "HyQMOM.jl"))

# Import the main.jl from src directory
include(joinpath(@__DIR__, "HyQMOM.jl", "src", "main.jl"))

