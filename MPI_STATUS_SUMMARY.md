# MPI Implementation Status

## ✅ What Works Perfectly

1. **Halo exchange mechanism**: Bitwise identical data transfer
2. **Unified main() interface**: Single entry point for serial/MPI
3. **Cartesian domain decomposition**: Correct subdomain assignment
4. **Serial vs MPI-1 match**: 0.000e+00 difference (BITWISE IDENTICAL)

## ❌ Current Issue

**MPI ranks produce different results from each other** when using 2+ ranks.

## 🔍 Root Cause Analysis

### The Tradeoff

We face a fundamental tension:

**Option A: Pass interior+halos to pas_HLL**
- ✓ pas_HLL has stencil access across processor boundaries
- ✓ MPI ranks are self-consistent
- ✗ Array size differs from serial (24 vs 20 for Np=20)
- ✗ Serial vs MPI-1 gives 0.074 error

**Option B: Pass interior-only to pas_HLL** (current)
- ✓ Array size matches serial exactly
- ✓ Serial vs MPI-1 is bitwise identical
- ✗ MPI ranks differ from each other (0.384 error)
- ✗ pas_HLL can't access neighbors across processor boundaries

### The Stencil Problem

`pas_HLL` computes flux at cell interfaces using:
```matlab
lleft(j) = min([vpmin(j), vpmin(j+1)]);  % Needs NEXT cell
Wstar(j) = ... M(j) ... M(j+1) ...       % Needs NEXT cell  
```

For a processor boundary at cell nx:
- Serial: M(nx+1) is cell nx+1 of the global domain
- MPI: M(nx+1) would be in the HALO (from neighbor rank)
- Current: We don't pass halos, so processor boundaries get physical BCs instead

## 💡 Possible Solutions

### Solution 1: Modify pas_HLL for MPI
Add boundary-type awareness to pas_HLL:
- Accept arrays with halos
- Know which boundaries are physical vs processor
- Apply BCs only at physical boundaries
- **Issue**: Already tried, still had 0.047 error

### Solution 2: Split pas_HLL calls
- Physical boundaries: pass Np cells, apply physical BCs
- Processor boundaries: pass Np+1 cells (interior + 1 halo), skip BCs
- **Complexity**: Different logic per boundary type

### Solution 3: Manual flux computation at processor boundaries
- Use pas_HLL for bulk interior
- Manually compute fluxes at processor boundaries using halo data
- **Best approach**: Most explicit control

### Solution 4: Increase halo to include stencil
- Current halo=2, but stencil only needs 1
- Revert halo to 1, ensure communication is correct
- **Issue**: Already tested, still had errors

## 📊 Test Results Summary

| Configuration | Difference | Status |
|---------------|-----------|---------|
| Serial (Np=20) vs Serial (Np=20) | 0.000e+00 | ✅ |
| Serial vs MPI-1 | 0.000e+00 | ✅ |
| MPI-1 vs MPI-2 | 3.844e-01 | ❌ |
| MPI-1 vs MPI-4 | 7.410e-01 | ❌ |

## 🎯 Recommendation

Implement **Solution 3**: Manual boundary flux computation
- This gives explicit control over processor vs physical boundaries
- Matches how finite volume methods typically handle domain decomposition
- Clear separation of concerns

