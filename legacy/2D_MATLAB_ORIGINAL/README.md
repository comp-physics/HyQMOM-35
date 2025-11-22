# 2D_MATLAB_ORIGINAL - Original 2D HyQMOM Code (Archive)

## Overview

This directory contains the **original, unorganized 2D HyQMOM MATLAB code** as it was developed during initial research. This code is preserved for historical reference and algorithm verification, but **should not be used for active development**.

## Purpose

This directory serves as:

1. **Historical archive** - Original research code as developed
2. **Algorithm reference** - Detailed mathematical implementations
3. **Verification baseline** - Cross-checking organized implementations
4. **Documentation** - Comments and notes from original development

## Important Notes

- **DO NOT use this for new development** - Use `../2D_MATLAB/` instead
- **Code is unorganized** - All files in flat directory structure
- **No formal test suite** - Some test files but not comprehensive
- **Inconsistent naming** - Function names vary in style and convention
- **Limited documentation** - Inline comments may be sparse or outdated

## Why Keep This?

While the code is unorganized, it contains:

- **Original algorithm implementations** that may have useful details
- **Intermediate development files** showing algorithm evolution
- **Test cases** that validate mathematical correctness
- **Comments and derivations** from original research

## What's Here?

This directory contains a flat collection of `.m` files including:

### Core Algorithm Files
- `hyqmom_3D.m` - Original HyQMOM algorithm
- `collision35.m` - BGK collision operator
- `delta2star3D.m`, `delta2star3D_permutation.m` - Delta-2 star operations
- `delta2starchol_L3.m` - Cholesky decomposition variant

### Moment Operations
- `InitializeM4_35.m` - Moment initialization
- `M4toC4_3D.m`, `C4toM4_3D.m` - Moment/cumulant conversions
- `C5toM5_3D.m` - 5th-order conversions
- `M2CS4_35.m` - Moment to C and S matrices
- `S4toC4_3D.m`, `S4toC4_3D_r.m` - S to C matrix conversions
- `Moments5_3D.m` - 5th-order moment calculations

### Realizability Checking
- `realizable_2D.m`, `realizable_3D.m` - General realizability
- `realizability_S2.m` - Second-order realizability
- `realizability_S111.m`, `realizability_S210.m`, `realizability_S211.m` - Higher-order checks
- `realizability_S310.m`, `realizability_S310_220.m`, `realizablity_S220.m` - Specific order checks
- `check2D.m` - 2D realizability verification
- `lower_bound_S220.m`, `bound_minor1.m` - Realizability bounds

### Flux and Numerics
- `Flux_closure35_and_realizable_3D.m` - Flux closure with realizability
- `flux_HLL.m` - HLL Riemann solver
- `pas_HLL.m` - HLL positive-advective scheme
- `closure_and_eigenvalues.m` - Closure and eigenvalue computation
- `eigenvalues6x_hyperbolic_3D.m`, `eigenvalues6y_hyperbolic_3D.m` - Eigenvalue calculations
- `jacobian6.m` - Jacobian matrix computation

### Utilities and Helpers
- `rootsR.m`, `rootsR_X_Y.m` - Root finding for quadrature
- `edge_corner_correction.m` - Boundary corrections

### Main Simulation and Visualization
- `main_2Dcrossing_3DHyQMOM35.m` - Main simulation script (2D crossing jets)
- `plot_final_time.m` - Final state visualization
- `plot3Dsym_C.m`, `plot3Dsym_mom.m`, `plot3Dsym_S.m` - 3D symmetric plotting
- `Cmoment3D_plots.m`, `Smoment3D_plots.m` - Moment field plots
- `contour3D_plots.m` - 3D contour plots
- `hyperbolic3D_plots.m` - Hyperbolic system visualization

### Test and Validation
- `test_symmetry_2D.m` - Symmetry verification

### Data Files
- `riemann_3D_hyqmom35_crossing_Np40_Kn1_Ma0.mat` - Saved simulation results

## Organized Alternative

**For active development, use the organized version in `../2D_MATLAB/`**, which contains:

- Proper `src/` directory structure
- Comprehensive `tests/` suite
- Consistent function interfaces
- MPI support with proper halo exchange
- Golden file generation
- Documentation and README

## Using This Code

If you need to reference something from this directory:

### 1. Check the algorithm details
```matlab
% Open the file to understand the mathematical implementation
edit hyqmom_3D.m
```

### 2. Compare with organized version
```matlab
% Check if organized version matches
diff('2D_MATLAB_ORIGINAL/hyqmom_3D.m', '2D_MATLAB/src/hyqmom_3D.m')
```

### 3. Extract useful functions
If you find something useful, port it to the organized codebase:
- Copy to `../2D_MATLAB/src/` or `../3D_MATLAB/src/`
- Update function name and interface for consistency
- Add to test suite
- Document properly

## Running Original Code (Not Recommended)

If you must run the original code:

```matlab
% Add all files to path
addpath('/path/to/2D_MATLAB_ORIGINAL');

% Run main simulation
main_2Dcrossing_3DHyQMOM35
```

**Warning**: The original code may:
- Have hard-coded paths
- Expect specific working directories
- Have undocumented dependencies
- Produce inconsistent results

## Key Differences from Organized Code

| Aspect | Original (this dir) | Organized (2D_MATLAB) |
|--------|---------------------|----------------------|
| Structure | Flat, all files in one directory | `src/`, `tests/`, proper hierarchy |
| Naming | Inconsistent | Standardized conventions |
| Testing | Minimal, ad-hoc | Comprehensive test suite |
| MPI | Limited or none | Full support with halo exchange |
| Documentation | Sparse comments | README + inline docs |
| Maintenance | Archived, read-only | Active development |

## Migration History

The organized code in `../2D_MATLAB/` was created by:

1. **Extracting core algorithms** from this directory
2. **Restructuring** into `src/` and `tests/`
3. **Standardizing** function interfaces
4. **Adding** MPI parallelization
5. **Creating** test suite and golden files
6. **Documenting** all components

## When to Reference This

Reference this original code when:

- Verifying mathematical correctness of algorithms
- Understanding original research intent
- Debugging subtle numerical differences
- Exploring algorithm variations tried during research
- Writing papers that cite original implementation

Do NOT use this for:

- New feature development
- Production simulations
- Teaching examples
- Collaborative projects

## Historical Context

This code represents the exploratory phase of HyQMOM development where:

- Algorithms were being actively developed and tested
- Various approaches were tried and compared
- Mathematical correctness was prioritized over code organization
- Files were created as needed without formal structure

The "mess" here is a natural part of research code evolution!

## Relationship to Other Implementations

- **2D_MATLAB** - Organized version of this code (use this instead)
- **3D_MATLAB** - Extended to 3D with similar organization
- **HyQMOM.jl** - Julia port with modern software practices
- **HyQMOM_2D_archive** - Archived Julia port of 2D code

## Citation

If referencing the original implementation:

```bibtex
@software{hyqmom_original,
  title = {HyQMOM: Original 2D Implementation},
  author = {Spencer H. Bryngelson and contributors},
  year = {2024},
  note = {Archived research code}
}
```

## Contact

For questions about:
- **Algorithm details**: Reference this directory, check inline comments
- **Using the code**: Use `../2D_MATLAB/` instead with its README
- **Historical context**: Contact repository maintainers

