# 2D_MATLAB - 2D HyQMOM MATLAB Implementation

## Overview

This directory contains the organized, production-ready 2D Hyperbolic Quadrature Method of Moments (HyQMOM) implementation in MATLAB. This is the cleaned-up version of the original 2D code, structured similarly to the 3D implementation with proper organization, testing, and MPI support.

## Purpose

The 2D implementation serves multiple purposes:

1. **Development platform** for testing algorithmic improvements before 3D implementation
2. **Faster testing** with reduced computational cost (2D vs 3D)
3. **Reference implementation** for 2D crossing jet problems
4. **Educational tool** for understanding HyQMOM with simpler geometry

## Differences from 2D_MATLAB_ORIGINAL

This directory (`2D_MATLAB/`) is the **organized, production version** with:

- Proper directory structure (`src/`, `tests/`)
- Consistent function naming and interfaces
- Comprehensive test suite
- MPI parallelization support
- Golden file generation for validation
- Documentation and comments

For the original, unorganized code, see `../2D_MATLAB_ORIGINAL/`.

## Quick Start

### Prerequisites

- MATLAB R2020b or later
- Parallel Computing Toolbox (optional, for MPI)
- MEX compiler configured (run `mex -setup`)

### Running a Simulation

```matlab
% Add paths
addpath(genpath('src'));

% Run basic simulation
main
```

### Building MEX Files

```matlab
% Build performance-critical MEX files
build_mex
```

## Features

- **2D moment-based solver** for Boltzmann-BGK equation
- **MPI parallelization** with 1D domain decomposition (x-direction)
- **HLL flux** with realizability preservation
- **Crossing jet test problem** as standard validation case
- **Golden file generation** for regression testing

## Project Structure

```
2D_MATLAB/
├── README.md                    # This file
├── main.m                       # Main simulation entry point
├── setup_paths.m                # Path configuration
├── simulation_plots.m           # Visualization
├── build_mex.m                  # MEX compilation
├── create_goldenfiles.m         # Golden file generation
├── src/                         # Source code
│   ├── apply_flux_update.m      # 2D flux application
│   ├── compute_halo_fluxes_and_wavespeeds.m
│   ├── halo_exchange_2d.m       # MPI halo exchange
│   ├── setup_mpi_cartesian_2d.m # 2D MPI setup
│   ├── hyqmom_3D.m              # Core algorithm (3D version)
│   ├── collision35.m            # Collision operator
│   ├── realizability.m          # Realizability checking
│   ├── autogen/                 # Auto-generated code
│   └── ...
├── tests/                       # Test suite
│   ├── run_all_tests.m          # Master test runner
│   ├── test_mpi_goldenfile.m    # MPI consistency
│   └── ...
└── goldenfiles/                 # 2D-specific golden files
    └── goldenfile_mpi_*ranks_*.mat
```

## Running Tests

### Full Test Suite

```matlab
cd tests
exit_code = run_all_tests();
```

### Individual Tests

```matlab
% MPI consistency test
test_mpi_goldenfile

% Local MPI test
test_mpi_local_only
```

### Creating Golden Files

```matlab
% CI mode: 1-2 ranks (fast)
create_goldenfiles('ci')

% Local mode: 4-8 ranks (comprehensive)
create_goldenfiles('local')
```

## MPI Parallelization

### Domain Decomposition (2D)

- **x-direction**: Divided among MPI ranks (1D decomposition)
- **y-direction**: Full data on all ranks (no decomposition)
- **Halo cells**: Ghost cells exchanged between neighbors

### Running with MPI

```matlab
% In main.m, configure:
mpi_enabled = true;
num_ranks = 2;

% Then run
main
```

## Test Problem: Crossing Jets

The standard 2D test case simulates two jets crossing perpendicular to each other:

- **Domain**: 2D spatial domain (x, y)
- **Initial condition**: Two Gaussian jet profiles
- **Physics**: Boltzmann-BGK collision operator
- **Parameters**: Configurable Knudsen and Mach numbers

## Relationship to Other Implementations

### vs. 2D_MATLAB_ORIGINAL
- This is the **organized version** with proper structure
- Original code in `../2D_MATLAB_ORIGINAL/` is kept for historical reference
- Use this version for active development

### vs. 3D_MATLAB
- Similar structure and organization
- Simpler geometry (2D vs 3D) for faster testing
- Some functions shared (e.g., `hyqmom_3D.m` works for both)

### vs. HyQMOM_2D_archive (Julia)
- `../HyQMOM_2D_archive/` contains archived Julia port of 2D code
- Current Julia work focuses on 3D in `../HyQMOM.jl/`

## Development Workflow

1. **Test new algorithms** in 2D (faster iteration)
2. **Validate** with golden files and test suite
3. **Port to 3D** once validated in 2D
4. **Cross-check** 2D and 3D implementations

## Performance Tips

### Quick Tests
```matlab
Np = 20;      % Small grid
tmax = 0.1;   % Short time
CFL = 0.7;
```

### Production Runs
```matlab
Np = 100;     % High resolution
tmax = 1.0;   % Full evolution
CFL = 0.5;    % Conservative stability
num_ranks = 4; % Parallel
```

## Troubleshooting

### Path Issues
```matlab
% Ensure paths are set
setup_paths();

% Or manually
addpath(genpath('src'));
```

### MPI Not Working
- Check Parallel Computing Toolbox: `ver('parallel')`
- Fall back to serial: `mpi_enabled = false`

### Golden File Mismatches
- Regenerate with current code: `create_goldenfiles('ci')`
- Check consistent parameters across runs

## Some Key Algorithms

### Moment System
Evolves 35 moments even in 2D (for consistency with 3D):
- Moment (0,0,0): Mass
- Moments (1,0,0), (0,1,0), (0,0,1): Momentum
- Moments (2,0,0), (1,1,0), ... : Higher-order moments

### Spatial Discretization
- Finite volume method on Cartesian grid
- HLL (Harten-Lax-van Leer) flux for hyperbolic terms
- Realizability-preserving limiter

### Time Integration
- Forward Euler with operator splitting
- Hyperbolic step + collision step
- CFL-limited time stepping

## Further Reading

- See `../3D_MATLAB/README.md` for 3D implementation details
- See `../2D_MATLAB_ORIGINAL/README.md` for original code history
- See `../HyQMOM.jl/README.md` for Julia implementation

## Contact

For questions:
- Check test suite for usage examples
- See main repository README
- Review `2D_MATLAB_ORIGINAL/` for algorithm background
