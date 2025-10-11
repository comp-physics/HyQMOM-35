# 2D MATLAB Code Archive

## Summary

The `2D_MATLAB/` directory contains the complete **original 2D MATLAB implementation** before the 3D physical space conversion.

## Quick Reference

- **Location**: `2D_MATLAB/`
- **Archived from**: Commit `2daa0cb` (October 10, 2024)
- **Last 2D commit before**: `f05d49f` ("matlab now 3d")
- **Total files**: 71 MATLAB files
- **Status**: Complete, self-contained, validated implementation

## What's Archived

| Component | Files | Description |
|-----------|-------|-------------|
| Main solver | `main.m` | 2D physical space driver |
| Simulation runner | `src/simulation_runner.m` | 2D time-stepping (3D moments: nx×ny×Nmom) |
| MPI utilities | `src/setup_mpi_cartesian_2d.m` | 2D domain decomposition |
| Halo exchange | `src/halo_exchange_2d.m` | 2D neighbor communication |
| Flux updates | `src/apply_flux_update.m` | X and Y direction updates |
| Tests | `tests/` | Complete 2D test suite (16 test files) |
| Plotting | `simulation_plots.m` | 2D visualization |

## Key Implementation Details

### Array Dimensions (2D version)
```matlab
M = zeros(nx+2*halo, ny+2*halo, Nmom);  % 3D array: 2D grid + moments
```

### Current 3D Version
```matlab
M = zeros(nx+2*halo, ny+2*halo, nz, Nmom);  % 4D array: 3D grid + moments
```

## Why Archive?

1. **Reference**: Validated 2D implementation for comparison
2. **Performance**: 2D is faster for certain studies
3. **Validation**: Cross-check 3D results with 2D (homogeneous-z cases)
4. **Documentation**: Track codebase evolution
5. **Research**: Some studies only need 2D physical space

## Usage Example

```matlab
cd 2D_MATLAB
addpath('src')
addpath('src/autogen')

% Run 2D simulation: 40×40 grid, 4 MPI ranks, time=0.1
results = main(40, 0.1, true, 4);

% Run tests
cd tests
run_all_tests
```

## Comparison: 2D vs 3D

| Feature | 2D (archived) | 3D (current) |
|---------|---------------|--------------|
| Physical grid | `Np×Np` | `Np×Np×Nz` |
| Moment array | 3D | 4D |
| MPI decomposition | 2D (x,y) | 2D (x,y only, full z) |
| Initial conditions | Jets in x-y plane | Jets extruded in z |
| Parameters | `Np, tmax` | `Np, tmax, Nz, homogeneous_z` |
| Test suite | 16 tests | 18 tests (adds 3D tests) |

## Files NOT Archived

- `original/` - Legacy unstructured code (unchanged)
- `goldenfiles/` - Generated data (can be regenerated)
- `HyQMOM.jl/` - Julia implementation (separate development)
- Git history - Available via `git show 2daa0cb:path/to/file`

## Notes

- This is a **static snapshot**, not actively maintained
- For active development, use the main 3D code
- To access any file from the original 2D commit: `git show 2daa0cb:path/to/file`
- Tests in this archive expect 2D golden files (not included)

---

*Archive created: October 11, 2025*
*Purpose: Preserve validated 2D MATLAB implementation*

