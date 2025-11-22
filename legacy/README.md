# Legacy Code Archive

This directory contains archived implementations and auxiliary code from the HyQMOM project. These files are preserved for historical reference, cross-validation, and algorithm verification.

## Directory Structure

### MATLAB Implementations

#### `3D_MATLAB/` - 3D MATLAB Reference Implementation
- **Purpose**: Production 3D MATLAB code used for cross-validation with Julia
- **Status**: Maintained for golden file generation and numerical reference
- **Usage**: See [3D_MATLAB/README.md](3D_MATLAB/README.md)
- **Key Features**:
  - MPI parallelization support
  - Golden file generation for Julia validation
  - Comprehensive testing suite
  - Serves as numerical reference for Julia port

#### `2D_MATLAB/` - Organized 2D MATLAB Implementation
- **Purpose**: Production 2D MATLAB code with proper structure
- **Status**: Maintained for algorithm development and validation
- **Usage**: See [2D_MATLAB/README.md](2D_MATLAB/README.md)
- **Key Features**:
  - Faster testing platform (2D vs 3D)
  - Algorithm development and validation
  - MPI support with 1D decomposition
  - Clean, organized codebase

#### `2D_MATLAB_ORIGINAL/` - Original 2D Research Code
- **Purpose**: Original unorganized 2D MATLAB research code
- **Status**: Archived, preserved for historical reference
- **Usage**: See [2D_MATLAB_ORIGINAL/README.md](2D_MATLAB_ORIGINAL/README.md)
- **Note**: Not for active development (use `2D_MATLAB/` instead)

### Julia Archives

#### `HyQMOM_2D_archive/` - Archived 2D Julia Port
- **Purpose**: Early 2D Julia implementation
- **Status**: Archived, superseded by current 3D HyQMOM.jl
- **Usage**: See [HyQMOM_2D_archive/README.md](HyQMOM_2D_archive/README.md)
- **Note**: Kept for historical reference and algorithm verification

### Simulation Outputs

#### `snapshots/`, `old-snaps/`, `new-snaps/` - Archived Simulation Data
- **Purpose**: Historical simulation outputs and visualization snapshots
- **Status**: Archived
- **Contents**: PNG images, JLD2 snapshot files, visualization data
- **Note**: For reference only; not part of active development

#### `for-rodney/` - Specific Simulation Results
- **Purpose**: Simulation results prepared for collaboration/sharing
- **Status**: Archived
- **Contents**: Selected snapshots and visualization data

#### `HyQMOM.jl-snapshots/` - Development Artifacts
- **Purpose**: Miscellaneous development files from the Julia package
- **Status**: Archived
- **Contents**: Experimental visualization code, temporary snapshots, development notes

## Active Development

**The current, actively maintained Julia implementation is at the repository root**, not in this `legacy/` directory.

For active development and production use:
- Source code: `../src/`
- Tests: `../test/`
- Examples: `../examples/`
- Documentation: `../docs/`

## Cross-Validation Workflow

The MATLAB implementations in this directory are still used for:

1. **Golden File Generation**: Creating reference data for Julia validation
   ```matlab
   cd legacy/3D_MATLAB
   create_goldenfiles('ci')  % Creates files in ../goldenfiles/
   ```

2. **Numerical Reference**: Verifying Julia results against MATLAB
   ```bash
   # From repository root
   julia --project=. test/test_golden_files.jl
   ```

3. **Algorithm Development**: Testing new methods in 2D before 3D implementation
   ```matlab
   cd legacy/2D_MATLAB
   setup_paths()
   main
   ```

## Maintenance Notes

- **MATLAB code**: Updated only for bug fixes or golden file regeneration
- **Julia archives**: Frozen; no updates expected
- **Simulation outputs**: Preserved as-is for reference
- **CI/CD**: Workflows in `../.github/workflows/` reference these paths

## Questions?

For questions about:
- **Legacy code**: Review the individual README files in each subdirectory
- **Active development**: See the main repository README at `../README.md`
- **Contributing**: Focus on the main Julia package, not legacy code

---

**Last Updated**: November 2025  
**Restructure**: Repository reorganized to promote HyQMOM.jl package to root

