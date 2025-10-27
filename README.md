# HyQMOM - Hyperbolic Quadrature Method of Moments

3D Hyperbolic Quadrature Method of Moments (HyQMOM) solver for the Boltzmann equation with BGK collision operator, featuring MPI parallelization and interactive visualization.

## Repository Structure

This repository contains multiple implementations and versions of the HyQMOM solver:

### **[HyQMOM.jl/](HyQMOM.jl/)** - Main Julia Implementation (Current/Active)
The production-ready 3D Julia implementation with:
- Interactive 3D visualization with GLMakie
- MPI parallelization with 2D domain decomposition
- Comprehensive test suite and examples

**Start here for new development and production runs.**

See [HyQMOM.jl/README.md](HyQMOM.jl/README.md) for detailed documentation.

#### 📚 Documentation
Build and view comprehensive documentation locally:
```bash
# Build documentation
./build_docs.sh        # Linux/macOS
# Serve documentation locally  
./serve_docs.sh        # Starts server at http://localhost:8000
```

### **[3D_MATLAB/](3D_MATLAB/)** - 3D MATLAB Implementation
Reference 3D MATLAB implementation used for cross-validation:
- MPI parallelization support
- Golden file generation for Julia validation
- Comprehensive testing suite
- Serves as numerical reference for Julia port

See [3D_MATLAB/README.md](3D_MATLAB/README.md) for usage instructions.

### **[2D_MATLAB/](2D_MATLAB/)** - 2D MATLAB Implementation (Organized)
Production 2D MATLAB code with proper structure:
- Faster testing platform (2D vs 3D)
- Algorithm development and validation
- MPI support with 1D decomposition
- Clean, organized codebase

See [2D_MATLAB/README.md](2D_MATLAB/README.md) for details.

### **[2D_MATLAB_ORIGINAL/](2D_MATLAB_ORIGINAL/)** - Original 2D Code (Archive)
Original unorganized 2D MATLAB research code:
- Historical reference and algorithm verification
- Preserved for mathematical correctness checking
- Not for active development (use 2D_MATLAB instead)

See [2D_MATLAB_ORIGINAL/README.md](2D_MATLAB_ORIGINAL/README.md) for context.

### **[HyQMOM_2D_archive/](HyQMOM_2D_archive/)** - Archived 2D Julia Port
Archived 2D Julia implementation:
- Superseded by current 3D HyQMOM.jl
- Kept for historical reference and algorithm verification
- Not actively maintained

See [HyQMOM_2D_archive/README.md](HyQMOM_2D_archive/README.md) for details.

## Quick Start

### Julia (Recommended)

```bash
cd HyQMOM.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run interactive 3D visualization
julia --project=. examples/run_3d_jets_timeseries.jl

# With MPI parallelization
mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60
```

### MATLAB (3D Reference)

```matlab
cd 3D_MATLAB
setup_paths()
main  % Run simulation
```

## Golden Files

The `goldenfiles/` directory at the repository root contains MATLAB-generated reference files (`.mat` format) used for cross-validation between MATLAB and Julia implementations.

### Generating Golden Files

```matlab
% From MATLAB (3D_MATLAB/)
cd 3D_MATLAB
create_goldenfiles('ci')  % Creates files in ../goldenfiles/
```

### Validating Julia Against MATLAB

```bash
# From Julia (HyQMOM.jl/)
cd HyQMOM.jl
julia --project=. test/test_golden_files.jl
```

The golden files ensure numerical agreement between implementations (typically to machine precision ~1e-15).

## Testing

### Julia Tests

```bash
cd HyQMOM.jl

# Full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# MPI tests
./test/run_mpi_tests.sh

# All tests including MPI
./test/run_tests.sh
```

### MATLAB Tests

```matlab
% 3D MATLAB tests
cd 3D_MATLAB/tests
run_all_tests()

% 2D MATLAB tests
cd 2D_MATLAB/tests
run_all_tests()
```

## CI/CD

The repository includes GitHub Actions workflows for:

- **Julia validation** against MATLAB golden files
- **MATLAB test suite** execution
- **MPI consistency** testing
- **Cross-implementation** verification

See `.github/workflows/` for configuration details.

## Features

### Physical Model
- **Boltzmann equation** with BGK collision operator
- **3D moment-based** kinetic solver (35 moments)
- **Hyperbolic quadrature** method for moment inversion
- **Realizability preservation** at all stages

### Numerical Methods
- **Finite volume** spatial discretization
- **HLL flux** (Harten-Lax-van Leer) with realizability
- **Forward Euler** time integration with operator splitting
- **CFL-limited** time stepping

### Parallelization
- **MPI domain decomposition**: 2D (x-y plane) for 3D code
- **Automatic halo exchange** for ghost cells
- **Scales** to multiple nodes on HPC systems

### Visualization (Julia)
- **Interactive 3D** isosurface rendering with GLMakie
- **Time-series animation** with play/pause controls
- **Multiple quantities**: density, velocity, pressure, temperature
- **Mouse controls**: rotate, zoom, pan

## Performance

### Quick Testing
```bash
# Julia: small grid, short time
julia --project=HyQMOM.jl HyQMOM.jl/examples/run_3d_jets_timeseries.jl --Nx 20 --Ny 20 --tmax 0.01
```

### Production Runs
```bash
# Julia: high resolution, MPI parallel
mpiexec -n 8 julia --project=HyQMOM.jl HyQMOM.jl/examples/run_3d_jets_timeseries.jl --Nx 100 --Ny 100 --tmax 0.5
```

## Documentation

Each subdirectory contains detailed README files:

| Directory | Description | README |
|-----------|-------------|--------|
| `HyQMOM.jl/` | Current Julia implementation | [README](HyQMOM.jl/README.md) |
| `3D_MATLAB/` | 3D MATLAB reference | [README](3D_MATLAB/README.md) |
| `2D_MATLAB/` | Organized 2D MATLAB | [README](2D_MATLAB/README.md) |
| `2D_MATLAB_ORIGINAL/` | Original 2D archive | [README](2D_MATLAB_ORIGINAL/README.md) |
| `HyQMOM_2D_archive/` | Archived 2D Julia | [README](HyQMOM_2D_archive/README.md) |

## Development Workflow

1. **Algorithm development**: Test in 2D MATLAB (faster iteration)
2. **Validation**: Generate golden files and verify
3. **3D implementation**: Port to 3D MATLAB or Julia
4. **Cross-validation**: Compare implementations using golden files
5. **Production**: Use Julia for performance and scalability

## System Requirements

### Julia
- Julia 1.9 or later
- MPI.jl (automatically configured)
- GLMakie (for visualization, optional in CI)

### MATLAB
- MATLAB R2020b or later
- Parallel Computing Toolbox (for MPI)
- MEX compiler configured

## HPC Usage

### Slurm Batch Scripts

The repository includes example SLURM batch scripts for HPC systems:

```bash
# Base template
sbatch hyqmom_base.sbatch

# Specific simulation example
sbatch hyqmom_Ma0p0_Kn1p0_t0p04.sbatch
```

Edit the `.sbatch` files to adjust:
- Number of nodes and tasks
- Grid resolution (`--Nx`, `--Ny`, `--Nz`)
- Simulation time (`--tmax`)
- Project account

## License

See [license.md](license.md) for licensing information.

## Contact

For questions, issues, or contributions:
- Open an issue on GitHub
- See implementation-specific READMEs for detailed usage
- Check test suites for API usage examples

---

## Quick Reference Card

### Julia Commands
```bash
# Install dependencies
cd HyQMOM.jl && julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run simulation with visualization
julia --project=HyQMOM.jl HyQMOM.jl/examples/run_3d_jets_timeseries.jl

# Run with MPI (4 ranks)
mpiexec -n 4 julia --project=HyQMOM.jl HyQMOM.jl/examples/run_3d_jets_timeseries.jl --Nx 60 --Ny 60

# Run tests
cd HyQMOM.jl && julia --project=. -e 'using Pkg; Pkg.test()'
```

### MATLAB Commands
```matlab
% 3D MATLAB
cd 3D_MATLAB
setup_paths()
main

% Generate golden files
create_goldenfiles('ci')

% Run tests
cd tests
run_all_tests()
```

### Common Parameters
```bash
--Nx 40           # Grid resolution in x direction
--Ny 40           # Grid resolution in y direction
--Nz 40           # Grid resolution in z direction
--tmax 0.2        # Simulation time
--Ma 1.0          # Mach number
--Kn 1.0          # Knudsen number
--CFL 0.7         # CFL number
```

---

**For detailed usage, see [HyQMOM.jl/README.md](HyQMOM.jl/README.md)** and the USERGUIDE.
