# MPI Golden File Testing

## Overview

Golden file testing for MPI implementation with 1, 2, 3, and 4 ranks.
All configurations use **20 grid points per rank per dimension** with **tmax = 0.1**.

## Grid Configurations

| Ranks | Process Grid | Grid Size | Points/Rank |
|-------|--------------|-----------|-------------|
| 1     | 1×1          | 20×20     | 20×20       |
| 2     | 2×1          | 40×40     | 20×20       |
| 3     | 3×1          | 60×60     | 20×20       |
| 4     | 2×2          | 40×40     | 20×20       |

## Quick Start

### 1. Create Golden Files

```matlab
create_goldenfiles
```

This will create:
- `goldenfiles/goldenfile_mpi_1ranks_Np20_tmax100.mat`
- `goldenfiles/goldenfile_mpi_2ranks_Np40_tmax100.mat`
- `goldenfiles/goldenfile_mpi_3ranks_Np60_tmax100.mat`
- `goldenfiles/goldenfile_mpi_4ranks_Np40_tmax100.mat`

### 2. Run Tests

```matlab
runtests('tests/test_mpi_goldenfile.m')
```

Or run individual tests:
```matlab
runtests('tests/test_mpi_goldenfile.m/test_mpi_1_rank')
runtests('tests/test_mpi_goldenfile.m/test_mpi_2_ranks')
runtests('tests/test_mpi_goldenfile.m/test_mpi_3_ranks')
runtests('tests/test_mpi_goldenfile.m/test_mpi_4_ranks')
```

## Test Details

### test_mpi_1_rank
- Grid: 20×20
- Tests MPI implementation with single rank
- Should match golden file within ~1e-6

### test_mpi_2_ranks
- Grid: 40×40 (2×1 decomposition)
- Tests 2-rank MPI decomposition
- Should match golden file within ~1e-6

### test_mpi_3_ranks
- Grid: 60×60 (3×1 decomposition)
- Tests 3-rank MPI decomposition
- Should match golden file within ~1e-6

### test_mpi_4_ranks
- Grid: 40×40 (2×2 decomposition)
- Tests 2D MPI decomposition (2×2)
- Should match golden file within ~1e-6

### test_mpi_consistency_across_ranks
- Compares results between different rank counts
- Informational test (different grids = different physics)

## Expected Results

All tests should pass with maximum differences ~1e-6, which is expected
due to floating-point non-associativity in parallel reductions.

## Files

- `create_goldenfiles.m` - Creates golden files for all rank configurations
- `tests/test_mpi_goldenfile.m` - Unit tests for MPI implementation
- `main_mpi.m` - MPI implementation with domain decomposition

## Requirements

- MATLAB R2023b or later
- Parallel Computing Toolbox
- Minimum grid size: 10 points per rank per dimension
