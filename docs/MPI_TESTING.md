# MPI Implementation Testing Guide

This document describes the testing strategy for the MPI domain decomposition implementation.

## Overview

The MPI implementation is validated through:
1. **Consistency tests** - MPI results match serial reference
2. **Golden file tests** - Regression testing across code changes
3. **Rank consistency tests** - Results consistent across 1, 2, and 4 ranks

## Test Files

### Core Test Files

| File | Purpose |
|------|---------|
| `tests/test_mpi_goldenfile.m` | Main MPI test suite with golden file validation |
| `main_mpi.m` | MPI-parallel version of solver |
| `run_mpi_goldenfile_creation.m` | Script to create MPI golden files |

### Golden Files

| File | Description |
|------|-------------|
| `goldenfiles/goldenfile_Np10_tmax100.mat` | Serial reference (from `main.m`) |
| `goldenfiles/goldenfile_mpi_1ranks_Np10_tmax100.mat` | MPI with 1 rank |
| `goldenfiles/goldenfile_mpi_2ranks_Np10_tmax100.mat` | MPI with 2 ranks |
| `goldenfiles/goldenfile_mpi_4ranks_Np10_tmax100.mat` | MPI with 4 ranks |

## Running Tests

### 1. Create Golden Files

First, generate the MPI golden files:

```matlab
run_mpi_goldenfile_creation
```

This will:
- Create golden files for 1, 2, and 4 ranks
- Verify symmetry preservation
- Save files to `goldenfiles/` directory

Expected output:
```
═══════════════════════════════════════════════════════════
Creating golden file for 1 rank(s)...
═══════════════════════════════════════════════════════════

Running MPI simulation with 1 rank(s)...
Simulation completed in XX.XX seconds
  Final time: 0.100000
  Time steps: XXX

✓ Golden file saved: goldenfiles/goldenfile_mpi_1ranks_Np10_tmax100.mat
...
```

### 2. Run Test Suite

Execute the full MPI test suite:

```matlab
% Run all tests
runtests('tests/test_mpi_goldenfile.m')

% Run specific test
runtests('tests/test_mpi_goldenfile.m', 'Name', 'test_mpi_1_rank_vs_serial')
```

### 3. Individual Tests

Run specific validation tests:

```matlab
% Test 1 rank vs serial
runtests('tests/test_mpi_goldenfile.m/test_mpi_1_rank_vs_serial')

% Test 2 ranks vs golden
runtests('tests/test_mpi_goldenfile.m/test_mpi_2_ranks_vs_golden')

% Test 4 ranks vs golden
runtests('tests/test_mpi_goldenfile.m/test_mpi_4_ranks_vs_golden')

% Test consistency across ranks
runtests('tests/test_mpi_goldenfile.m/test_mpi_consistency_across_ranks')
```

## Test Descriptions

### Test 1: MPI 1 Rank vs Serial
**Purpose:** Verify that MPI with 1 rank produces identical results to serial code

**What it tests:**
- MPI infrastructure doesn't alter results
- Halo exchange works correctly even with no neighbors
- Boundary conditions applied correctly

**Tolerance:** 10⁻¹² (strict)

**Expected result:**
```
=== TEST: MPI 1 Rank vs Serial ===
Running MPI simulation with 1 rank...
  Completed: Np=10, tmax=0.100, ranks=1, steps=XXX
Comparing: 1-rank MPI vs Serial
  Final time diff: X.XXXe-XX (tolerance: 1.000e-09) ✓
  Time steps diff: 0 (max: 5) ✓
  M: max=X.XXe-XX, mean=X.XXe-XX ✓
  C: max=X.XXe-XX, mean=X.XXe-XX ✓
  S: max=X.XXe-XX, mean=X.XXe-XX ✓
  M5: max=X.XXe-XX, mean=X.XXe-XX ✓
  C5: max=X.XXe-XX, mean=X.XXe-XX ✓
  S5: max=X.XXe-XX, mean=X.XXe-XX ✓
  PASS: 1-rank MPI vs Serial
```

### Test 2: MPI 2 Ranks vs Golden
**Purpose:** Regression test for 2-rank MPI implementation

**What it tests:**
- 2-rank decomposition remains stable across code changes
- Halo exchange between 2 subdomains
- Load balancing with 2 ranks

**Tolerance:** 10⁻¹²

### Test 3: MPI 4 Ranks vs Golden
**Purpose:** Regression test for 4-rank MPI implementation

**What it tests:**
- 4-rank decomposition (2×2 grid) remains stable
- Corner halo exchanges work correctly
- Load balancing with 4 ranks

**Tolerance:** 10⁻¹²

### Test 4: Consistency Across Ranks
**Purpose:** Verify results are independent of rank count

**What it tests:**
- 1, 2, and 4 ranks produce identical results
- Domain decomposition doesn't introduce errors
- Communication patterns preserve accuracy

**Comparisons:**
- 1 rank vs 2 ranks
- 1 rank vs 4 ranks
- 2 ranks vs 4 ranks

**Tolerance:** 10⁻¹²

## Understanding Test Results

### Successful Test Output

```
Running test_mpi_goldenfile

=== TEST: MPI 1 Rank vs Serial ===
...
  PASS: 1-rank MPI vs Serial

=== TEST: MPI 2 Ranks vs Golden File ===
...
  PASS: 2-rank MPI vs Golden

=== TEST: MPI 4 Ranks vs Golden File ===
...
  PASS: 4-rank MPI vs Golden

=== TEST: MPI Consistency Across Ranks ===
...
  PASS: 4-rank vs 1-rank

Done test_mpi_goldenfile
__________

All 4 tests passed.
```

### Error Indicators

**High errors (> 10⁻¹²):**
- Likely issue with halo exchange
- Check boundary condition application
- Verify domain decomposition logic

**Dimension mismatches:**
- Problem with gather operation
- Check subdomain indexing
- Verify block partitioning

**Time step differences (> 5):**
- CFL condition might differ slightly
- Check global reduction in time step calculation
- Verify wave speed calculations

## Tolerance Guidelines

| Quantity | Tolerance | Reason |
|----------|-----------|--------|
| Moments (M, C, S, M5, C5, S5) | 10⁻¹² | Strict for physics |
| Final time | 10⁻⁹ | Allow small CFL variations |
| Time steps | ±5 steps | Floating point in dt calculation |

## Troubleshooting

### Test Failure: "Golden file not found"

**Solution:**
```matlab
run_mpi_goldenfile_creation
```

### Test Failure: High errors in moments

**Possible causes:**
1. **Halo exchange issue**
   - Check `halo_exchange_2d.m`
   - Verify neighbor topology
   - Test with 1 rank first

2. **Boundary condition problem**
   - Check `apply_physical_bc_2d.m`
   - Verify global boundary detection
   - Compare with serial BC application

3. **Gather operation error**
   - Check `main_mpi.m` gather logic
   - Verify subdomain indexing
   - Print intermediate values

**Debug steps:**
```matlab
% Enable verbose output in main_mpi
% Add breakpoints in halo_exchange_2d
% Run with 1 rank first, then 2, then 4
```

### Test Failure: Dimension mismatch

**Check:**
1. Block partitioning logic
2. Halo width accounting
3. Interior vs total array sizes

### Performance Issues

**If tests are slow:**
1. Reduce grid size for development:
   ```matlab
   % In test file, temporarily use smaller grid
   Np = 6;  % instead of 10
   ```

2. Use fewer iterations:
   ```matlab
   tmax = 0.05;  % instead of 0.1
   ```

3. Test with 1 rank only for quick validation

## CI/CD Integration

### GitHub Actions

Add to `.github/workflows/test.yml`:

```yaml
- name: Run MPI Golden File Tests
  run: |
    matlab -batch "addpath('tests'); results = runtests('test_mpi_goldenfile'); assertSuccess(results)"
```

### Pre-commit Hook

Create `.git/hooks/pre-commit`:

```bash
#!/bin/bash
matlab -batch "runtests('tests/test_mpi_goldenfile')" || {
    echo "MPI tests failed. Commit aborted."
    exit 1
}
```

## Verification Checklist

Before pushing MPI code changes:

- [ ] All 4 MPI tests pass
- [ ] Errors < 10⁻¹² for all moment fields
- [ ] 1-rank matches serial exactly
- [ ] 2-rank and 4-rank consistent
- [ ] Golden files up to date
- [ ] Symmetry preserved (< 10⁻¹²)
- [ ] No dimension mismatches
- [ ] Time steps within ±5 of expected

## Golden File Structure

Each golden file contains:

```matlab
golden_data = 
  struct with fields:
    parameters: [struct]  % Np, tmax, num_workers, final_time, time_steps, etc.
    grid: [struct]        % x, y, xm, ym, dx, dy
    moments: [struct]     % M, C, S, M5, C5, S5
    metadata: [struct]    % creation_date, matlab_version, num_ranks, etc.
    filename: [char]
```

## Performance Metrics

Expected test times (on typical workstation):

| Test | Ranks | Grid | Time | 
|------|-------|------|------|
| 1 rank vs serial | 1 | 10×10 | ~30s |
| 2 ranks vs golden | 2 | 10×10 | ~30s |
| 4 ranks vs golden | 4 | 10×10 | ~30s |
| Consistency | 1,2,4 | 10×10 | ~90s |

**Total test suite time:** ~3 minutes

## Advanced Testing

### Custom Test Parameters

Modify test parameters for specific scenarios:

```matlab
% In test file, add custom configuration
function test_mpi_custom(testCase)
    % Test with different grid size
    mpi_data = run_mpi_simulation(16, 0.1, 4);
    % Custom validation...
end
```

### Stress Testing

Test with extreme configurations:

```matlab
% Large grid
main_mpi(64, 0.1, false, 4)

% Many ranks
main_mpi(32, 0.05, false, 16)

% Long integration
main_mpi(10, 1.0, false, 4)
```

## References

- Main test file: `tests/test_mpi_goldenfile.m`
- MPI solver: `main_mpi.m`
- Serial reference: `main.m`
- MPI utilities: `src/setup_mpi_cartesian_2d.m`, `src/halo_exchange_2d.m`

---

*Last updated: October 1, 2025*

