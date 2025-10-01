# CI Testing with MPI Implementation

## Issue: MPI Tests Fail in GitHub Actions

### Problem
The GitHub Actions CI was failing because MPI tests require the **Parallel Computing Toolbox (PCT)**, which is not included in the standard MATLAB CI license.

**Error:**
```
Undefined function 'gcp' for input arguments of type 'char'.
```

Functions requiring PCT:
- `gcp()` - Get current parallel pool
- `parpool()` - Create parallel pool
- `spmd` - Single Program Multiple Data
- `labSendReceive()` - Inter-worker communication
- `labBroadcast()` - Broadcast to all workers

### Solution
Make MPI tests **conditional** - they skip gracefully when PCT is unavailable.

## Implementation

### Detection Logic
```matlab
% In setupOnce:
testCase.TestData.has_pct = license('test', 'Distrib_Computing_Toolbox') && ...
                             ~isempty(ver('parallel'));
```

This checks:
1. **License available:** `license('test', 'Distrib_Computing_Toolbox')`
2. **Toolbox installed:** `~isempty(ver('parallel'))`

### Skip Logic (Each Test)
```matlab
function test_mpi_1_rank_vs_serial(testCase)
    % Skip if Parallel Computing Toolbox not available
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 1 Rank vs Serial ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    % ... rest of test
end
```

**Key points:**
- Use `assumeFail()` not `error()` - marks test as "skipped" not "failed"
- Clear messaging for why test was skipped
- Early return prevents execution of PCT-dependent code

## CI Behavior

### Before Fix
```
Running test_mpi_goldenfile
Error: Undefined function 'gcp'

Failure Summary:
  4 Failed, 0 Incomplete
  test_mpi_goldenfile/test_mpi_1_rank_vs_serial    X  Errored
  test_mpi_goldenfile/test_mpi_2_ranks_vs_golden   X  Errored
  test_mpi_goldenfile/test_mpi_4_ranks_vs_golden   X  Errored
  test_mpi_goldenfile/test_mpi_consistency_across_ranks X Errored

❌ CI FAILS
```

### After Fix
```
Running test_mpi_goldenfile

=== TEST: MPI 1 Rank vs Serial ===
SKIPPED: Parallel Computing Toolbox not available

=== TEST: MPI 2 Ranks vs Golden File ===
SKIPPED: Parallel Computing Toolbox not available

=== TEST: MPI 4 Ranks vs Golden File ===
SKIPPED: Parallel Computing Toolbox not available

=== TEST: MPI Consistency Across Ranks ===
SKIPPED: Parallel Computing Toolbox not available

Summary:
  0 Passed, 0 Failed, 4 Incomplete (assumed)

Running test_validate_goldenfile
=== RIGOROUS GOLDEN FILE VALIDATION ===
... (serial test runs normally)
  
Summary:
  1 Passed, 0 Failed

✅ CI PASSES
```

## Testing Environments

| Environment | PCT Available | MPI Tests | Serial Tests | CI Status |
|-------------|---------------|-----------|--------------|-----------|
| **GitHub Actions CI** | ❌ No | Skipped | ✅ Run | ✅ Pass |
| **Local Dev (with PCT)** | ✅ Yes | ✅ Run | ✅ Run | ✅ Pass |
| **Local Dev (no PCT)** | ❌ No | Skipped | ✅ Run | ✅ Pass |

## Local Development

### With Parallel Computing Toolbox
All tests run normally:
```matlab
>> runtests('tests/test_mpi_goldenfile.m')

Running test_mpi_goldenfile
=== TEST: MPI 1 Rank vs Serial ===
... (full test execution)
  PASS: 1-rank MPI vs Serial

=== TEST: MPI 2 Ranks vs Golden File ===
  PASS: 2-rank MPI vs Golden

=== TEST: MPI 4 Ranks vs Golden File ===
  PASS: 4-rank MPI vs Golden

=== TEST: MPI Consistency Across Ranks ===
  PASS: 4-rank vs 1-rank

All 4 tests passed.
```

### Without Parallel Computing Toolbox
Tests skip gracefully:
```matlab
>> runtests('tests/test_mpi_goldenfile.m')

Warning: MPI tests require Parallel Computing Toolbox - tests will be skipped

Running test_mpi_goldenfile
=== TEST: MPI 1 Rank vs Serial ===
SKIPPED: Parallel Computing Toolbox not available

... (all 4 tests skipped)

0 Passed, 0 Failed, 4 Incomplete.
```

## Files Modified

**Commit:** `b0e2986`

**Changed:**
- `tests/test_mpi_goldenfile.m` (+44 lines)
  - Added PCT detection in `setupOnce()`
  - Added skip logic to all 4 test functions
  - Used `assumeFail()` for graceful skipping

## Verification

### Check PCT Availability
```matlab
% Check license
has_license = license('test', 'Distrib_Computing_Toolbox')

% Check installation
has_toolbox = ~isempty(ver('parallel'))

% Both required
has_pct = has_license && has_toolbox
```

### Manual Test
```matlab
% Should work if PCT available
pool = gcp('nocreate')

% Should work regardless
runtests('tests/test_validate_goldenfile.m')

% Behavior depends on PCT
runtests('tests/test_mpi_goldenfile.m')
```

## Best Practices for CI

### ✅ Do This
- **Separate test files** for serial and parallel tests
- **Conditional execution** based on toolbox availability
- **Use `assumeFail()`** for skipped tests (not `error()`)
- **Clear messaging** about why tests skipped
- **Document requirements** in test header

### ❌ Don't Do This
- **Don't** assume all toolboxes available in CI
- **Don't** use `error()` for missing dependencies
- **Don't** mix required and optional tests without checks
- **Don't** hide skip messages

## Future Improvements

### Option 1: Separate Test Files
```
tests/
  test_validate_goldenfile.m      (serial, always runs)
  test_mpi_goldenfile.m           (parallel, conditional)
```

### Option 2: Test Tags
```matlab
% Mark tests with tags
tests = functiontests(localfunctions);
tests(1).Tags = {'MPI', 'RequiresPCT'};

% Run only tests without PCT requirement
runtests('tests', 'Tag', '~RequiresPCT')
```

### Option 3: CI-Specific Workflow
```yaml
# .github/workflows/test.yml
- name: Run Serial Tests Only
  if: ${{ !contains(github.event.head_commit.message, '[test-mpi]') }}
  run: matlab -batch "runtests('tests/test_validate_goldenfile.m')"
```

## Related Documentation

- **MPI Testing Guide:** `docs/MPI_TESTING.md`
- **MPI Testing Summary:** `docs/MPI_TESTING_SUMMARY.md`
- **MPI README:** `docs/MPI_README.md`

## References

- MATLAB Unit Testing: `matlab.unittest.TestCase`
- Assumptions: `assumeFail()`, `assumeTrue()`, `assumeEqual()`
- License checking: `license('test', 'toolbox_name')`
- Version checking: `ver('toolbox')`

---

**Status:** ✅ CI now passes with MPI tests gracefully skipping when PCT unavailable

**Commit:** `b0e2986` - Fix CI: Skip MPI tests when Parallel Computing Toolbox unavailable

