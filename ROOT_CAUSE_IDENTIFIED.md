# ROOT CAUSE IDENTIFIED

## The Bug

**File:** `src/halo_exchange_2d.m` and wave speed exchange in `main_mpi.m`
**Issue:** Halo exchange was hardcoded for halo=1, sending only a single slice

### Original (buggy) code:
```matlab
labSend(A(h+1:h+1, h+1:h+ny, :), left_neighbor);  % WRONG: sends 1 column regardless of h
```

### Fixed code:
```matlab
labSend(A(h+1:h+h, h+1:h+ny, :), left_neighbor);  % CORRECT: sends h columns
```

## Test Results

| Halo Width | 1v2 ranks | 1v3 ranks | Status |
|------------|-----------|-----------|--------|
| halo=1 (buggy) | 0.073 error | ❌ Failed | BUG |
| halo=1 (fixed) | 0.109 error | ✅ 0.0 error | PARTIAL |
| halo=2 (fixed) | ✅ 0.0 error | (not tested) | PERFECT |

## Why halo=1 still fails for 2 ranks

Looking at pas_HLL stencil:
- Interior update at j: uses F(j), F(j+1), M(j), M(j+1)
- With halo=1: Only 1 neighbor cell available
- With 2 decomposition: The interface is at a specific location where stencil may need 2 cells

For 3 ranks (1x3 grid), each rank has same size, decomposition is symmetric.
For 2 ranks (1x2 grid), asymmetry may cause stencil to need wider halo at specific locations.

**OR** the issue is lines 18-19 of pas_HLL that overwrite F boundaries!

## Conclusion

**Required halo width: 2 cells** (not 1)

This is likely because:
1. pas_HLL's stencil actually needs 2-cell access in some configurations
2. OR pas_HLL's boundary overwrites (lines 18-19) need special handling

**Fix:** Use halo=2 by default
