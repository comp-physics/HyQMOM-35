# 3D_MATLAB - 3D HyQMOM MATLAB Implementation

## Overview

This directory contains the 3D Hyperbolic Quadrature Method of Moments (HyQMOM) implementation in MATLAB for solving the Boltzmann equation with the BGK collision operator. This is the production MATLAB implementation with MPI parallelization support.

## Features

- **3D moment-based kinetic solver** for Boltzmann-BGK equation
- **MPI parallelization** with domain decomposition
- **Cartesian grid decomposition** in x-y plane
- **HLL flux** computation with realizability preservation
- **Golden file generation** for cross-validation with Julia implementation

## Quick Start

### Prerequisites

- MATLAB R2020b or later
- Parallel Computing Toolbox (for MPI support)
- MEX compiler configured (run `mex -setup` in MATLAB)

### Running a Simulation

```matlab
% Add paths
setup_paths();

% Run basic simulation
main
```

### Building MEX Files

Some performance-critical functions use MEX (compiled C code):

```matlab
% Build all MEX files
build_mex
```

## MPI Parallelization

### Running with MPI

The implementation supports MPI parallelization through MATLAB's Parallel Computing Toolbox:

```matlab
% Configure MPI settings in main.m
mpi_enabled = true;
num_ranks = 4;

% Run simulation
main
```

### Domain Decomposition

- **x-y plane**: Divided among MPI ranks using 2D Cartesian decomposition
- **z-direction**: Full data replicated on all ranks
- **Halo exchange**: Automatic ghost cell communication between neighboring ranks

## Golden File Generation

Golden files are used to validate the Julia implementation against MATLAB results. These files are stored in the root `goldenfiles/` directory (shared with Julia tests).

### Creating Golden Files

```matlab
% CI mode: 1 and 2 ranks (fast, for automated testing)
create_goldenfiles('ci')

% Local mode: 4 and 8 ranks (comprehensive testing)
create_goldenfiles('local')

% All configurations
create_goldenfiles('all')
```

Golden files are saved to `../../goldenfiles/` (repository root) with naming convention:
- `goldenfile_mpi_1ranks_Np20_tmax100.mat`
- `goldenfile_mpi_2ranks_Np20_tmax100.mat`
- etc.

## Testing

### Running Test Suite

```matlab
cd tests
run_all_tests()
```

### Test Categories

- **Unit tests**: Individual function validation (`test_*.m`)
- **MPI tests**: Multi-rank consistency verification (`test_mpi_*.m`)
- **Integration tests**: Full simulation validation
- **Golden file tests**: Verification against reference data

### Test Files

- `run_all_tests.m` - Master test runner
- `test_mpi_goldenfile.m` - MPI consistency testing
- `test_mpi_local_only.m` - Local MPI verification
- `create_test_goldenfiles.m` - Test-specific golden file generation

## Project Structure

```
3D_MATLAB/
├── README.md                    # This file
├── main.m                       # Main simulation entry point
├── setup_paths.m                # Path configuration
├── simulation_plots.m           # Visualization functions
├── build_mex.m                  # MEX compilation script
├── create_goldenfiles.m         # Golden file generation
├── src/                         # Source code
│   ├── apply_flux_update_3d.m
│   ├── compute_halo_fluxes_and_wavespeeds_3d.m
│   ├── halo_exchange_3d.m
│   ├── setup_mpi_cartesian_3d.m
│   ├── hyqmom_3D.m              # Core HyQMOM algorithm
│   ├── collision35.m            # BGK collision operator
│   ├── realizability.m          # Realizability checks
│   ├── autogen/                 # Auto-generated symbolic code
│   └── ...                      # Additional modules
└── tests/                       # Test suite
    ├── run_all_tests.m
    ├── test_mpi_goldenfile.m
    └── ...
```

## Key Algorithms

### Moment Evolution

The code evolves 35 moments of the distribution function:
- Mass, momentum (3), pressure tensor (6), heat flux (3), higher-order moments (22)

### HyQMOM Algorithm

Core quadrature method in `src/hyqmom_3D.m`:
1. Moment inversion to find quadrature points and weights
2. Realizability checking and correction
3. Flux computation at cell interfaces
4. Conservative moment update

### Collision Operator

BGK collision operator in `src/collision35.m`:
- Relaxation toward local Maxwellian
- Preserves mass, momentum, and energy

## Relationship to Julia Implementation

This MATLAB code serves as the reference implementation for the Julia port (now at repository root). The two implementations:

- **Share golden files**: Both use `../../goldenfiles/*.mat` for cross-validation
- **Use identical algorithms**: Julia port is a direct translation
- **Maintain numerical agreement**: Typically to machine precision (~1e-15)

To verify Julia against MATLAB:
```matlab
% Generate golden files from MATLAB
create_goldenfiles('ci')

% Then run Julia tests (from repository root)
% cd ../HyQMOM.jl
% julia --project=. test/test_golden_files.jl
```

## Performance Tips

### For Quick Testing
- Use small grid sizes (Np=20)
- Short time spans (tmax=0.1)
- Single rank

### For Production Runs
- Increase grid resolution (Np=100+)
- Use multiple MPI ranks
- Compile MEX files for critical functions

### Memory Optimization
- Monitor memory usage with `memory` command
- Reduce grid size if needed
- Use fewer MPI ranks on limited-memory systems

## Troubleshooting

### MEX Compilation Issues
```matlab
% Reconfigure MEX compiler
mex -setup

% Try default compiler
mex -setup C
```

### MPI Not Available
- Ensure Parallel Computing Toolbox is installed: `ver('parallel')`
- Run in serial mode: Set `mpi_enabled = false` in `main.m`

### Path Issues
- Always run `setup_paths()` before using functions
- Or add to MATLAB startup: `addpath(genpath('src'))`

### Golden File Mismatches
- Ensure same parameters (Np, tmax, etc.)
- Check MATLAB version compatibility
- Verify MEX files are properly compiled

## Citation

If you use this code in research, please cite:

```bibtex
@software{hyqmom_matlab,
  title = {HyQMOM: 3D Hyperbolic Quadrature Method of Moments},
  author = {Spencer H. Bryngelson and contributors},
  year = {2024},
  note = {MATLAB Implementation}
}
```

## Contact

For questions or issues:
- See main repository README
- Check Julia implementation in `../HyQMOM.jl/`
- Review test suite for usage examples

