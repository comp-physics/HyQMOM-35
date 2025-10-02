# MPI Domain Boundary Issue

## Problem
Results differ by ~0.07 per domain boundary between different MPI decompositions.
- 1 rank vs 2 ranks: max difference = 0.073
- 1 rank vs 4 ranks: max difference = 0.146 (exactly 2×)
- Each decomposition is self-consistent and reproducible

## Root Cause
The `pas_HLL` function applies boundary conditions at array endpoints:
```matlab
F(1,:) = F(2,:);     % Line 18
F(Np,:) = F(Np-1,:); % Line 19
```

In parallel:
- These endpoints are processor boundaries with valid halo data from neighbors
- But `pas_HLL` overwrites F(1) and F(Np) with BC values
- This causes the first/last interior cells to use incorrect flux values

## Why Fixing It Is Complex
1. Simply not applying BC makes it worse (0.09 instead of 0.07)
2. Conditionally applying BC based on neighbor flags also makes it worse
3. The interaction between Wstar BC (needed) and F BC (problematic) is subtle
4. The HLL flux computation couples M, F, and wave speeds in complex ways

## Impact
- Error scales linearly with number of domain boundaries
- ~0.3% relative error per boundary for this test case
- Acceptable for many applications but not bitwise reproducible

## Recommendations
1. For production: use consistent decomposition and validate with convergence studies
2. For exact reproducibility: use single rank or document tolerance
3. Long-term: redesign pas_HLL to handle processor boundaries explicitly

## Files Affected
- `src/pas_HLL.m`: BC application at lines 16-19
- `main_mpi.m`: Calls pas_HLL with exchanged halo data
