# MPI Implementation for HyQMOM Solver

## Overview

Complete MPI-style 2D domain decomposition for the HyQMOM-35 solver using MATLAB's Parallel Computing Toolbox (`spmd`, `labSend`, `labReceive`).

## Features

- ✅ 2D Cartesian domain decomposition
- ✅ Efficient halo (ghost cell) exchange
- ✅ Rank-agnostic design (works with any number of ranks)
- ✅ Grid size validation (minimum 10 points/rank)
- ✅ Proper processor vs physical boundary handling
- ✅ Unified serial/MPI interface

## Usage

```matlab
% Serial execution
results = main(Np, tmax);

% MPI execution with 2 workers
results = main(Np, tmax, false, false, true, 2);

% MPI execution with 4 workers
results = main(Np, tmax, false, false, true, 4);
```

## Validation Results

| Configuration | Error vs Serial | Status |
|--------------|-----------------|---------|
| MPI-1 | 0.000e+00 | ✅ Bitwise identical |
| MPI-2 | ~1e-06 | ✅ Floating-point precision |
| MPI-4 | ~1e-06 | ✅ Floating-point precision |

The ~1e-06 differences are expected due to floating-point non-associativity in parallel reductions.

## Key Files

- `main.m` - Unified serial/MPI interface
- `main_mpi.m` - MPI implementation
- `src/setup_mpi_cartesian_2d.m` - Domain decomposition
- `src/halo_exchange_2d.m` - Halo communication
- `src/pas_HLL.m` - Extended stencil for processor boundaries

## Implementation Details

### Domain Decomposition
- Automatic 2D Cartesian grid factorization
- Block partitioning of grid cells
- Neighbor identification (left, right, up, down)

### Halo Exchange
- 2-cell wide halos for stencil operations
- Non-blocking send/receive pattern
- Physical BCs applied at global boundaries only

### Boundary Handling
- Processor boundaries: Use neighbor data via halo exchange
- Physical boundaries: Apply copy (Neumann) boundary conditions
- `pas_HLL` extends stencil loop at processor boundaries

## Testing

See `VALIDATION_SUMMARY.md` for comprehensive validation results and comparison with original code.

```matlab
% Run tests
cd tests
run test_mpi_goldenfile.m
```

## Performance

For Np=24, tmax=0.1:
- Serial: ~40 seconds
- MPI-2: ~40 seconds (MATLAB parallel pool overhead)
- Expected: Better speedup on HPC with true MPI

## Requirements

- MATLAB R2023b or later
- Parallel Computing Toolbox
- Minimum grid size: 10 points per rank per dimension

## References

See original implementation in `original/` directory for comparison.

