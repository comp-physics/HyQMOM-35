# MPI Implementation - Debugging Status

## Current Achievement
✅ **Serial vs MPI-1: BITWISE IDENTICAL (0.000e+00)**

## Current Problem
❌ **MPI-2 and MPI-4 differ from Serial/MPI-1 by ~0.047**
- Error appears from **FIRST timestep** (~5.8e-7)
- Error at M(12,12,2) - right at the processor boundary!
- Error grows to ~0.047 by t=0.1

## What We've Tried

### Attempt 1: Pass interior-only to pas_HLL
- ✅ Result: Serial = MPI-1 (bitwise identical)
- ❌ Problem: MPI ranks differ (processor boundaries get physical BCs)

### Attempt 2: Compute wave speeds in halos
- Added vpxmin_ext, vpymin_ext with wave speeds for halo cells
- ❌ Still ~0.047 error

### Attempt 3: Recompute Fx/Fy in halos  
- Compute Fx, Fy, wave speeds in all halo cells after M exchange
- ❌ Still ~0.047 error

### Attempt 4: Conditional BC application in pas_HLL
- Pass apply_bc_left/right flags to skip BCs at processor boundaries
- ❌ Still ~0.047 error

## Root Cause Analysis

The fundamental issue: **`pas_HLL` stencil needs special handling at processor boundaries**

When we pass full array (interior+halos) to `pas_HLL`:
- At j=halo+1 (first interior cell on right rank), stencil needs M(j) and M(j+1)
- M(j+1) is at j=halo+2, which is interior data
- But pas_HLL also looks at M(1) and M(Np), which are halos
- The flux computation uses these halo values, but something is still wrong

## Next Steps

**Option A: Manual flux at boundaries**
Instead of passing halos to pas_HLL, keep interior-only but manually compute
the flux at the processor boundary using halo data.

**Option B: Debug pas_HLL stencil logic**
Figure out exactly why pas_HLL with halos+BC flags still gives wrong answer.

**Option C: Different HLL approach**
Maybe the issue is in how flux_HLL uses the data. Check if we need to 
handle the boundary flux computation differently.

