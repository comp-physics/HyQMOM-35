# HyQMOM_2D_archive - Archived 2D Julia Implementation

## Overview

This directory contains an **archived Julia implementation** of the 2D HyQMOM solver. This code represents an earlier Julia port that has been superseded by the current 3D implementation in `../HyQMOM.jl/`.

## Status: Archived

⚠️ **This code is archived and not actively maintained.**

For active Julia development, see: **`../HyQMOM.jl/`**

## Purpose

This archive serves as:

1. **Historical record** of 2D Julia port development
2. **Reference implementation** for 2D-specific features
3. **Comparison baseline** for algorithm validation
4. **Learning resource** for understanding Julia port evolution

## Why Archived?

Development effort has shifted to the 3D Julia implementation (`../HyQMOM.jl/`) because:

- **3D is more general** - Includes 2D as a special case
- **Better structure** - Improved package design and organization
- **Modern features** - Interactive visualization, better MPI support
- **Active maintenance** - Regular updates and improvements
- **Production ready** - Comprehensive testing and documentation

## What's Here?

This archive contains a Julia package structure:

```
HyQMOM_2D_archive/
├── Project.toml              # Julia package dependencies
├── Manifest.toml             # Dependency lock file
├── main.jl                   # Main simulation entry point
├── compare_goldenfiles.jl    # MATLAB comparison utilities
├── create_golden_for_comparison.jl  # Golden file generation
├── src/                      # Source code
│   ├── HyQMOM.jl            # Main module
│   ├── main.jl              # Simulation runner
│   ├── simulation_runner.jl # Core simulation loop
│   ├── moments/             # Moment operations
│   ├── realizability/       # Realizability checks
│   ├── numerics/            # Flux and collision operators
│   ├── mpi/                 # MPI parallelization
│   ├── utils/               # Utilities
│   ├── visualization/       # Plotting functions
│   └── autogen/             # Auto-generated symbolic code
└── test/                    # Test suite
    ├── test_mpi.jl
    ├── test_golden_files.jl
    └── ...
```

## Key Features (Historical)

When this was active, it provided:

- **2D Julia implementation** of HyQMOM algorithm
- **MPI parallelization** with 1D domain decomposition
- **Golden file validation** against MATLAB reference
- **Crossing jets test** problem
- **Auto-generated code** from symbolic computation

## Relationship to Current Implementation

### Evolution Path

1. **2D MATLAB** (`2D_MATLAB/`) - Original MATLAB implementation
2. **2D Julia** (this archive) - Initial Julia port
3. **3D Julia** (`HyQMOM.jl/`) - Current active development

### What Changed

From this 2D implementation to current 3D:

| Feature | 2D (Archived) | 3D (Current) |
|---------|---------------|--------------|
| Dimensions | 2D (x, y) | 3D (x, y, z) |
| Domain Decomposition | 1D (x-direction) | 2D (x-y plane) |
| Visualization | Basic PyPlot | Interactive GLMakie |
| Examples | Single main.jl | Multiple example scripts |
| Package Structure | Basic | Production-ready |
| Documentation | Limited | Comprehensive README |

### Code Similarities

Many algorithms are shared:
- Moment operations (InitializeM4_35, hyqmom_3D, etc.)
- Realizability checking
- HLL flux computation
- BGK collision operator

The 3D version extended these rather than replacing them.

## Using This Archive

### ⚠️ Not Recommended for Active Use

Use `../HyQMOM.jl/` instead. But if you need to reference this:

### Installation (Historical)

```bash
cd HyQMOM_2D_archive
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Running (Historical)

```bash
# Serial
julia --project=. main.jl

# MPI parallel
mpiexec -n 4 julia --project=. main.jl
```

### Comparing Against MATLAB

```bash
# Generate comparison with MATLAB golden files
julia --project=. compare_goldenfiles.jl
```

## Why Keep This Archive?

Despite being superseded, this code is valuable for:

### 1. Algorithm Verification
Cross-checking that 3D implementation generalizes correctly:
```julia
# Compare 2D slice of 3D results against 2D implementation
```

### 2. Performance Baseline
Historical performance metrics:
- 2D may be faster for specific 2D problems
- Useful for benchmarking optimization improvements

### 3. Educational Value
Simpler code structure for learning:
- Fewer dimensions easier to understand
- Good starting point for Julia beginners

### 4. Git History
Preserved in git for:
- Understanding design decisions
- Tracking algorithm evolution
- Referencing past bug fixes

## Migration Notes

If porting code from this archive to `../HyQMOM.jl/`:

### 1. Check if already exists
Many functions were already ported - check first:
```bash
grep -r "function_name" ../HyQMOM.jl/src/
```

### 2. Update for 3D
Add z-dimension handling:
```julia
# 2D version
M[ix, iy, :] = ...

# 3D version
M[ix, iy, iz, :] = ...
```

### 3. Update MPI calls
Domain decomposition changed from 1D to 2D:
```julia
# 2D: 1D decomposition
setup_mpi_cartesian_1d(comm, Np)

# 3D: 2D decomposition  
setup_mpi_cartesian_2d(comm, Np, Np)
```

### 4. Update tests
Add test to `../HyQMOM.jl/test/`:
```julia
@testset "Feature from 2D archive" begin
    # ...
end
```

## Comparison with Current HyQMOM.jl

To see what changed:

```bash
# Compare directory structures
diff -r HyQMOM_2D_archive/src/ HyQMOM.jl/src/

# Compare specific files
diff HyQMOM_2D_archive/src/HyQMOM.jl HyQMOM.jl/src/HyQMOM.jl
```

## Dependencies (Historical)

Original dependencies (see `Project.toml`):
- Julia 1.9+
- MPI.jl
- MAT.jl (MATLAB file I/O)
- PyPlot.jl (visualization)
- LinearAlgebra, Statistics (stdlib)

Current `HyQMOM.jl` uses:
- GLMakie (instead of PyPlot) for interactive 3D visualization
- Enhanced MPI capabilities
- Additional dependencies for 3D features

## Testing (Historical)

```bash
# Run full test suite
cd HyQMOM_2D_archive
julia --project=. -e 'using Pkg; Pkg.test()'

# Run specific tests
julia --project=. test/test_mpi.jl
```

**Note**: Tests may fail if dependencies are outdated.

## Known Issues

As archived code:

- ❌ May not work with latest Julia versions
- ❌ Dependencies may be outdated
- ❌ No active bug fixes
- ❌ Limited documentation
- ❌ Visualization may break on newer systems

## When to Use This vs. HyQMOM.jl

| Use Case | Use This Archive | Use HyQMOM.jl |
|----------|------------------|---------------|
| New development | ❌ | ✅ |
| Production runs | ❌ | ✅ |
| 2D-only problems | Maybe | ✅ (3D works for 2D) |
| Learning Julia port | ✅ | ✅ |
| Algorithm verification | ✅ | ✅ |
| Historical reference | ✅ | ❌ |
| Paper citations | Check date | ✅ (current) |

## Citation

If referencing this archived implementation:

```bibtex
@software{hyqmom_2d_julia_archive,
  title = {HyQMOM 2D: Archived Julia Implementation},
  author = {Spencer H. Bryngelson and contributors},
  year = {2024},
  note = {Archived - see HyQMOM.jl for current version},
  url = {https://github.com/[repo]/HyQMOM_2D_archive}
}
```

For current work, cite `HyQMOM.jl` instead.

## Further Reading

- **Current Julia implementation**: See `../HyQMOM.jl/README.md`
- **2D MATLAB reference**: See `../2D_MATLAB/README.md`
- **3D MATLAB reference**: See `../3D_MATLAB/README.md`
- **Original algorithms**: See `../2D_MATLAB_ORIGINAL/README.md`

## Contact

For questions about:
- **This archived code**: Check git history, compare with current version
- **Current Julia development**: See `../HyQMOM.jl/` and its README
- **Algorithm details**: See MATLAB implementations or current HyQMOM.jl

---

**Remember**: For active development, use `../HyQMOM.jl/` 🚀

