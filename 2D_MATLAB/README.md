# 2D MATLAB Code Archive

This directory contains the **original 2D MATLAB implementation** of the HyQMOM-35 solver, archived from commit `2daa0cb` (before the 3D physical space conversion).

## What's Here

- **Original 2D code**: Physical space is 2D (x, y), moment space is still 3D (u, v, w)
- **MPI parallelization**: 2D domain decomposition in x-y directions
- **Complete test suite**: All tests for the 2D implementation
- **Golden files**: Pre-computed reference solutions for validation

## Key Differences from Current 3D Code

| Feature | 2D Code (this archive) | 3D Code (main) |
|---------|------------------------|----------------|
| Physical space | 2D (x, y) | 3D (x, y, z) |
| Moment space | 3D (u, v, w) | 3D (u, v, w) |
| MPI decomposition | 2D (x, y) | 2D (x, y), no z decomposition |
| Initial conditions | Jets in x-y plane | Jets extruded in z-direction |
| Grid parameters | `Np` (x, y) | `Np` (x, y), `Nz` (z) |

## Usage

To run the 2D code:

```matlab
cd 2D_MATLAB
addpath('src')
addpath('src/autogen')
results = main(40, 0.1, true, 4);  % 40x40 grid, 4 workers
```

## Files Archived

- `main.m` - Main simulation driver (2D version)
- `simulation_plots.m` - Visualization for 2D results
- `create_goldenfiles.m` - Generate reference solutions
- `build_mex.m` - Build MEX files for performance
- `setup_paths.m` - Path configuration
- `src/` - Source code for 2D solver
- `tests/` - Test suite for 2D implementation

## Git Reference

- **Archived from**: Commit `2daa0cb` ("clean up")
- **Date**: Just before 3D conversion
- **Branch**: `julia3`

## Notes

This is a **static archive** for reference. For active development, use the main 3D code. This archive ensures the validated 2D implementation remains available for:

- Comparison and validation
- 2D-specific research
- Understanding the evolution of the codebase
- Benchmarking 2D vs 3D performance

---

*Created: October 11, 2025*
*Purpose: Preserve working 2D MATLAB implementation before 3D conversion*

