# Golden File Testing

## Overview

This directory contains rigorous regression tests that validate simulation results against a golden reference file.

## Test Configuration

**Golden File Test**: `test_validate_goldenfile.m`
- **Parameters**: Np=10, tmax=0.1
- **Timesteps**: 6 checkpoints
- **Tolerance**: 1e-10 (machine epsilon)
- **Runtime**: ~50-60 seconds

This configuration provides rigorous validation with longer time integration compared to smaller test cases.

## Golden File

**Location**: `../goldenfiles/goldenfile_Np10_tmax100.mat`

**Contents**:
- Moment arrays: M, C, S, M5, C5, S5
- Grid data: xm, ym
- Parameters: Np, tmax, final_time, time_steps
- Metadata: creation_date, MATLAB_version, git_branch

## Usage

### Running Tests

```matlab
% From tests/ directory
runtests('test_validate_goldenfile')

% Or from project root
cd tests
runtests
```

### Regenerating Golden File

If you need to regenerate the golden file (e.g., after intentional algorithm changes):

```matlab
% From project root
run('run_goldenfile_creation.m')
```

**⚠️ Warning**: Only regenerate the golden file after carefully verifying that:
1. Changes to results are intentional and correct
2. All MaxDiff values are zero (perfect symmetry)
3. You've documented the reason for regeneration

## Test Results

All tests should **PASS** with:
- `max diff: 0.00e+00` for all moment arrays
- Final time matches within 1e-6
- Time steps match within ±5 steps

## Validation Criteria

The test validates:
1. ✅ Final time convergence
2. ✅ Number of time steps
3. ✅ All moment arrays (M, C, S, M5, C5, S5) match exactly
4. ✅ Grid dimensions match
5. ✅ No NaN or Inf values

## Continuous Integration

This test is designed to work with `matlab-actions/run-tests@v2` for GitHub Actions CI/CD.

