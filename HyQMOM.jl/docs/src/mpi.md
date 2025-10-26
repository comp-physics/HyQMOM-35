# MPI & Parallelization

HyQMOM.jl is designed from the ground up for parallel execution using MPI (Message Passing Interface). This guide covers everything you need to know about running parallel simulations.

## Overview

HyQMOM.jl uses domain decomposition parallelization optimized for the 35-moment kinetic system:

- **xy-plane decomposition**: The computational domain is divided among MPI ranks in the x-y plane
- **z-direction replication**: Each rank holds the complete z-dimension data  
- **Moment-aware halo exchange**: All 35 moments synchronized between neighboring ranks
- **Realizability-preserving communication**: Ensures moment realizability across processor boundaries
- **Rank 0 visualization**: Interactive visualization runs only on the master rank (rank 0)
- **Load balancing**: Automatic distribution of grid points for optimal performance

## Quick Start

### Basic Parallel Execution

```bash
# Serial execution (1 process)
julia --project=. examples/run_3d_jets_timeseries.jl --Np 40

# Parallel execution (4 processes)
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100

# High-resolution parallel run
mpiexec -n 8 julia --project=. examples/run_3d_jets_timeseries.jl --Np 120 --Nz 60
```

### Automatic Detection

All HyQMOM.jl examples automatically detect whether they're running in serial or parallel mode:

- **Serial mode**: `MPI.Comm_size() == 1`
- **Parallel mode**: `MPI.Comm_size() > 1`

No code changes are needed to switch between serial and parallel execution.

## Domain Decomposition

### Spatial Decomposition Strategy

The computational domain is divided as follows:

```
Original domain: Nx × Ny × Nz
With P MPI ranks: Each rank gets (Nx/√P) × (Ny/√P) × Nz

Example with 4 ranks (2×2 decomposition):
Rank 0: [1:Nx/2,   1:Ny/2,   1:Nz]
Rank 1: [Nx/2+1:Nx, 1:Ny/2,   1:Nz]  
Rank 2: [1:Nx/2,   Ny/2+1:Ny, 1:Nz]
Rank 3: [Nx/2+1:Nx, Ny/2+1:Ny, 1:Nz]
```

### Ghost Cells and Halo Exchange

Each rank maintains ghost cells (halo regions) that store data from neighboring ranks:

- **Purpose**: Enable finite difference stencils across rank boundaries
- **Width**: Typically 2-3 cells depending on the numerical scheme
- **Synchronization**: Automatic exchange after each time step

### Load Balancing

For optimal performance, choose the number of MPI ranks such that:

1. **Perfect square decomposition**: Use ranks = 1, 4, 9, 16, 25, ... for best load balancing
2. **Sufficient work per rank**: Each rank should have at least 10×10 grid points
3. **Memory constraints**: More ranks = less memory per rank

## Performance Guidelines

### Choosing the Number of Ranks

**Good choices:**
```bash
# 1 rank: Serial execution
julia --project=. examples/run_3d_jets_timeseries.jl --Np 40

# 4 ranks: 2×2 decomposition
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 80

# 9 ranks: 3×3 decomposition  
mpiexec -n 9 julia --project=. examples/run_3d_jets_timeseries.jl --Np 120

# 16 ranks: 4×4 decomposition
mpiexec -n 16 julia --project=. examples/run_3d_jets_timeseries.jl --Np 160
```

**Avoid:**
```bash
# Non-square numbers (poor load balancing)
mpiexec -n 6 julia ...   # 6 doesn't factor nicely
mpiexec -n 10 julia ...  # 10 doesn't factor nicely

# Too many ranks for resolution (communication overhead)
mpiexec -n 16 julia ... --Np 40  # Only 2.5×2.5 per rank
```

### Scaling Guidelines

| Resolution | Recommended Ranks | Points per Rank | Use Case |
|------------|------------------|-----------------|----------|
| 40×40×40   | 1-4             | 40×40 to 20×20  | Development |
| 80×80×80   | 4-9             | 20×20 to 13×13  | Production |
| 120×120×120| 9-16            | 13×13 to 10×10  | High-res |
| 200×200×200| 16-25           | 12×12 to 10×10  | HPC |

### Memory Scaling

Memory usage per rank decreases approximately as:
```
Memory per rank ≈ Total memory / Number of ranks
```

Use more ranks to fit larger problems in available memory:

```bash
# High memory usage (single rank)
julia ... --Np 200  # ~200³ points on 1 rank

# Distributed memory (16 ranks)  
mpiexec -n 16 julia ... --Np 200  # ~50³ points per rank
```

## MPI Configuration

### Installation and Setup

**Ubuntu/Debian:**
```bash
sudo apt-get install mpich libmpich-dev
# or
sudo apt-get install openmpi-bin libopenmpi-dev
```

**macOS (Homebrew):**
```bash
brew install mpich
# or  
brew install open-mpi
```

**Verify installation:**
```bash
mpiexec --version
julia -e 'using MPI; println(MPI.versioninfo())'
```

### Julia MPI Configuration

HyQMOM.jl uses MPI.jl, which should automatically detect your system MPI:

```julia
using MPI
MPI.versioninfo()  # Check detected MPI implementation
```

If you encounter issues, you may need to rebuild MPI.jl:
```bash
julia -e 'using Pkg; Pkg.build("MPI")'
```

## Visualization in Parallel

### Automatic Rank 0 Visualization

Interactive visualization automatically runs only on rank 0:

```julia
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    # Launch GLMakie viewer
    interactive_3d_timeseries_streaming(filename, grid, params)
end
```

**Benefits:**
- No display requirements on compute nodes
- Reduced memory usage on worker ranks
- Centralized output file management

### Headless Parallel Execution

For HPC systems without display capabilities:

```bash
# Method 1: Environment variable
export HYQMOM_SKIP_PLOTTING=true
mpiexec -n 16 julia --project=. examples/run_3d_jets_timeseries.jl

# Method 2: Command-line flag
mpiexec -n 16 julia --project=. examples/run_3d_jets_timeseries.jl --no-viz true

# Method 3: CI detection (automatic)
export CI=true
mpiexec -n 16 julia --project=. examples/run_3d_jets_timeseries.jl
```

## Performance Scaling

### Memory and Communication Scaling

Each processor stores:
- **35 moments** × **local grid points** × **8 bytes/float64**
- **Halo cells**: Additional boundary data for neighbor communication
- **Temporary arrays**: For flux computation and realizability checking

**Memory estimate**: ~280 bytes × (Nx/√P) × (Ny/√P) × Nz per processor, where P is the number of processors.

**Communication Patterns:**

*Halo Exchange Frequency:*
- Every time step for spatial flux computation
- After realizability correction to maintain consistency
- During visualization data collection (rank 0 only)

*Communication Volume:*
- **Boundary faces**: 35 moments × 2 × (Nx/√P + Ny/√P) × Nz per time step
- **Corner exchanges**: 35 moments × 4 corners for diagonal neighbors

### Scaling Efficiency

**Strong scaling** (fixed problem size, more processors):
- **Ideal range**: 4-64 processors for typical problems
- **Communication overhead**: Becomes significant with >100 processors
- **Sweet spot**: 10,000-100,000 grid points per processor

**Weak scaling** (proportional problem size increase):
- **Linear scaling**: Up to 1000+ processors for large problems
- **Memory bandwidth**: Often the limiting factor on modern HPC systems

## HPC and Batch Systems

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=hyqmom
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=02:00:00
#SBATCH --mem=32G

module load julia/1.10
module load openmpi/4.1

export HYQMOM_SKIP_PLOTTING=true
export JULIA_NUM_THREADS=1

mpiexec -n 16 julia --project=. examples/run_3d_jets_timeseries.jl \
  --Np 200 --Nz 100 --tmax 1.0 --snapshot-interval 5
```

### PBS Example

```bash
#!/bin/bash
#PBS -N hyqmom_job
#PBS -l nodes=4:ppn=4
#PBS -l walltime=02:00:00
#PBS -l mem=64gb

cd $PBS_O_WORKDIR
module load julia openmpi

export HYQMOM_SKIP_PLOTTING=true

mpiexec -n 16 julia --project=. examples/run_3d_jets_timeseries.jl \
  --Np 160 --tmax 0.5 --no-viz true
```

## Debugging Parallel Runs

### Common MPI Issues

**MPI not found:**
```bash
# Check MPI installation
which mpiexec
mpiexec --version

# Rebuild MPI.jl if needed
julia -e 'using Pkg; Pkg.build("MPI")'
```

**Rank mismatch errors:**
```bash
# Ensure consistent Julia environments across nodes
julia --project=. -e 'using Pkg; Pkg.status()'

# Check MPI configuration
julia -e 'using MPI; MPI.Init(); println("Rank: ", MPI.Comm_rank(MPI.COMM_WORLD))'
```

**Communication timeouts:**
```bash
# Try fewer ranks
mpiexec -n 2 julia ...

# Check network connectivity between nodes
mpiexec -n 4 hostname
```

### Debugging Tools

**Print rank information:**
```julia
using MPI
MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)
size = MPI.Comm_size(MPI.COMM_WORLD)
println("Rank $rank of $size on $(gethostname())")
```

**Profile parallel performance:**
```bash
# Use Julia's built-in profiler on each rank
julia --project=. -e '
using Profile
@profile include("examples/run_3d_jets_timeseries.jl")
Profile.print()
'
```

## Advanced Topics

### Custom Domain Decomposition

For non-standard domain shapes or specialized hardware, you can modify the decomposition in `src/mpi/domain_decomposition.jl`.

### Communication Optimization

- **Overlap computation and communication**: Use non-blocking MPI calls
- **Minimize synchronization points**: Reduce `MPI.Barrier()` calls
- **Optimize data layout**: Ensure contiguous memory access patterns

### Hybrid Parallelization

Combine MPI with Julia threading:

```bash
export JULIA_NUM_THREADS=4
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl
# Total parallelism: 4 MPI ranks × 4 threads = 16 parallel workers
```

## Performance Monitoring

### Timing Analysis

```julia
# Built into examples - check output for timing information
julia --project=. examples/run_3d_jets_timeseries.jl --Np 80

# Look for output like:
# "Simulation completed in 45.2 seconds"
# "Time per step: 0.12 seconds"
# "Communication overhead: 5.3%"
```

### Memory Usage

```bash
# Monitor memory during execution
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl &
watch -n 1 'ps aux | grep julia'
```

### Scaling Studies

Test parallel efficiency:

```bash
# Weak scaling (constant work per rank)
mpiexec -n 1 julia ... --Np 40   # 40×40 per rank
mpiexec -n 4 julia ... --Np 80   # 40×40 per rank  
mpiexec -n 9 julia ... --Np 120  # 40×40 per rank

# Strong scaling (constant total work)
mpiexec -n 1 julia ... --Np 120  # 120×120 total
mpiexec -n 4 julia ... --Np 120  # 60×60 per rank
mpiexec -n 9 julia ... --Np 120  # 40×40 per rank
```
