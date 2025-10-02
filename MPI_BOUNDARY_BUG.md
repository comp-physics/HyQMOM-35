# MPI Boundary Bug - Identified Issue

## Problem
`main_mpi.m` produces different results than serial `main.m` even with 1 MPI rank.

## Symptoms
- Max difference: ~0.074 at boundary cells (e.g., cell (1,10,5))
- Interior cells match perfectly
- MPI ranks are bitwise identical to each other (0.000e+00 difference)
- Serial vs MPI-1: 0.074 difference at boundaries

## Root Cause (Suspected)
The MPI code applies physical boundary conditions differently than serial, even when there are no processor boundaries (1 rank case).

## Location
Likely in how `pas_HLL` applies BCs at lines 18-19, or how physical boundaries are handled in the halo exchange when there are no processor neighbors.

## Next Steps
1. Compare how `main.m` applies BCs vs `main_mpi.m`
2. Check if `apply_physical_bc_2d` is being called correctly
3. Verify that with 1 rank, the MPI code should behave identically to serial

## Workaround
For now, the unified `main()` interface works. Use:
- `main(..., false, false, false)` for serial (correct)
- `main(..., false, false, true, N)` for MPI with N ranks (boundaries differ slightly)

MPI consistency is perfect (bitwise identical across ranks), but serial vs MPI has ~0.3% relative error at boundaries.
