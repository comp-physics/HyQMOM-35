# MPI Golden File Testing - Quick Reference

## 🎯 Quick Start

### 1. Create Golden Files
```matlab
run_mpi_goldenfile_creation
```

### 2. Run All Tests
```matlab
runtests('tests/test_mpi_goldenfile.m')
```

### 3. Expected Output
```
Running test_mpi_goldenfile
Test 1/4: test_mpi_1_rank_vs_serial ✓
Test 2/4: test_mpi_2_ranks_vs_golden ✓
Test 3/4: test_mpi_4_ranks_vs_golden ✓
Test 4/4: test_mpi_consistency_across_ranks ✓

All 4 tests passed.
```

---

## 📋 What Gets Tested

| Test | Description | Validates |
|------|-------------|-----------|
| **1-rank vs Serial** | MPI with 1 rank vs `main.m` | MPI infrastructure is correct |
| **2-rank vs Golden** | MPI 2 ranks vs stored golden file | 2-rank decomposition stable |
| **4-rank vs Golden** | MPI 4 ranks vs stored golden file | 4-rank decomposition stable |
| **Rank Consistency** | Compare 1, 2, 4 ranks to each other | Results independent of ranks |

---

## 📁 Files Created

### Core Files
- `main_mpi.m` - MPI-parallel solver (300+ lines)
- `tests/test_mpi_goldenfile.m` - Test suite (290+ lines)
- `run_mpi_goldenfile_creation.m` - Golden file creator (130+ lines)
- `docs/MPI_TESTING.md` - Complete testing guide (600+ lines)

### Golden Files (created by script)
- `goldenfiles/goldenfile_mpi_1ranks_Np10_tmax100.mat`
- `goldenfiles/goldenfile_mpi_2ranks_Np10_tmax100.mat`
- `goldenfiles/goldenfile_mpi_4ranks_Np10_tmax100.mat`

---

## ✅ Success Criteria

All tests must:
- Pass with errors < **10⁻¹²** for moment fields
- Match time steps within **±5 steps**
- Preserve symmetry < **10⁻¹²**
- Complete in < **3 minutes** total

---

## 🔧 Troubleshooting

### Golden files missing
```matlab
run_mpi_goldenfile_creation  % Creates them automatically
```

### High errors (> 10⁻¹²)
1. Check halo exchange: `src/halo_exchange_2d.m`
2. Check boundary conditions: `src/apply_physical_bc_2d.m`
3. Run with 1 rank first to isolate issue

### Tests slow
```matlab
% Use smaller grid for testing
main_mpi(6, 0.05, false, 2)  % Instead of (10, 0.1, false, 4)
```

---

## 🚀 CI/CD Integration

Add to `.github/workflows/test.yml`:

```yaml
- name: Create MPI Golden Files
  run: matlab -batch "run_mpi_goldenfile_creation"

- name: Run MPI Tests
  run: matlab -batch "results = runtests('tests/test_mpi_goldenfile'); assertSuccess(results)"
```

---

## 📊 Test Coverage

| Component | Tested | How |
|-----------|--------|-----|
| Domain decomposition | ✅ | 1, 2, 4 rank tests |
| Halo exchange | ✅ | Consistency checks |
| Boundary conditions | ✅ | 1-rank vs serial |
| Load balancing | ✅ | Different rank counts |
| Gather operation | ✅ | All tests |
| Time stepping | ✅ | CFL consistency |
| Physics correctness | ✅ | vs serial golden file |

---

## 📈 Performance Expectations

| Configuration | Expected Time |
|---------------|---------------|
| 1 rank, 10×10 | ~30 seconds |
| 2 ranks, 10×10 | ~30 seconds |
| 4 ranks, 10×10 | ~30 seconds |
| **Full suite** | **~3 minutes** |

---

## 🎓 Usage Examples

### Run Specific Test
```matlab
runtests('tests/test_mpi_goldenfile.m', 'Name', 'test_mpi_1_rank_vs_serial')
```

### Run with Verbosity
```matlab
runtests('tests/test_mpi_goldenfile.m', 'OutputDetail', 'Verbose')
```

### Custom Parameters
```matlab
% Edit main_mpi.m defaults or:
results = main_mpi(16, 0.1, false, 4);  % 16×16 grid, 4 ranks
```

---

## 🔍 What Each Test Verifies

### Test 1: MPI 1-Rank vs Serial
- **Purpose:** Prove MPI gives same answer as serial
- **Critical for:** Validating MPI infrastructure
- **Checks:** All moments match to machine precision

### Test 2: MPI 2-Ranks vs Golden
- **Purpose:** Regression test for 2-rank decomposition
- **Critical for:** Detecting changes that break 2-rank case
- **Checks:** Results stable across code updates

### Test 3: MPI 4-Ranks vs Golden
- **Purpose:** Regression test for 4-rank decomposition (2×2 grid)
- **Critical for:** Detecting issues with corner exchanges
- **Checks:** Results stable across code updates

### Test 4: Rank Consistency
- **Purpose:** Verify answer doesn't depend on rank count
- **Critical for:** Proving domain decomposition is correct
- **Checks:** 1, 2, and 4 ranks all give same answer

---

## 📝 Commit Information

**Branch:** `mpi`  
**Commit:** `4eefe53`  
**Files Added:** 4  
**Lines Added:** 1,063  

**Test Status:** ✅ All tests passing (after golden file creation)

---

For complete documentation, see `docs/MPI_TESTING.md`

