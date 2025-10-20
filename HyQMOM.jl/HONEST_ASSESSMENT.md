
# Honest Assessment: 3D Unsplit Implementation Status

## Current State (October 20, 2025)

After extensive debugging and testing, here's the honest truth about the 3D unsplit implementation:

### ❌ **The 3D unsplit method is NOT working correctly**

## What We Discovered

### 1. **Original Implementation (commit 7442d44)**
- Had incorrect `flux_HLL` function calls
- Crashed with `MethodError: no method matching flux_HLL(...)`
- Never actually ran successfully

### 2. **After My "Fixes" (commits a61edb2 → 0078dc6)**
- Fixed `flux_HLL` call → prevented crashes
- Changed flux divergence indexing → **broke mass conservation**
- Split vs unsplit differ by ~4.5%
- Mass loss of ~0.85% over 20 timesteps

### 3. **Attempted Revert**
- Reverting flux indexing → crashes again (negative sqrt in collision)
- Original indexing was unstable
- My indexing breaks conservation
- **Neither version works correctly!**

## Root Problems

### Problem 1: Flux Indexing Confusion
The code mixes two different indexing conventions:
- **Flux storage**: `Fx[i]` = flux at face between cells i and i+1
- **Divergence**: Needs `(F_right - F_left) / dx` for cell i
- **Current state**: Inconsistent between X, Y, Z directions

### Problem 2: Z-Direction Has No Halo
- X and Y have halo cells for boundary conditions
- Z direction does NOT have halo
- This breaks the flux indexing pattern
- Boundary treatment is ad-hoc (`dFz = zeros(Nmom)` for boundaries)

### Problem 3: Never Validated Against MATLAB
- No golden file tests for unsplit method
- Dimensional splitting matches MATLAB perfectly
- Unsplit was never compared to known-good results
- We don't know what "correct" looks like for unsplit!

## What Actually Works

✅ **Dimensional Splitting Method**
- Matches MATLAB to machine precision  
- Mass conserved perfectly
- Stable and production-ready
- This is the validated, working code

## What We Thought We Had

We believed:
- ❌ Unsplit matched MATLAB (never tested)
- ❌ Split vs unsplit should agree (they don't, by ~4.5%)
- ❌ Rotational invariance of ~2% was acceptable (but method doesn't conserve mass!)

## The Hard Truth

**The 3D unsplit method needs to be implemented from scratch with:**

1. **Correct flux indexing from first principles**
   - Clear documentation of index conventions
   - Consistent treatment of all directions
   - Proper boundary conditions

2. **Validation against known results**
   - Golden file tests
   - Mass conservation verification
   - Agreement with dimensional splitting for non-rotated cases

3. **Z-direction halo implementation**
   - OR different boundary handling strategy
   - Must maintain conservation

## Recommendations

### Option 1: Start Fresh (Recommended)
1. Implement unsplit from scratch with clear design
2. Document flux indexing conventions
3. Add Z-halo support
4. Validate each step against dimensional splitting
5. Only then test rotational invariance

### Option 2: Accept Current State
1. Use dimensional splitting (working, validated)
2. Document that unsplit is experimental/broken
3. Accept ~20% rotational asymmetry from splitting
4. This is still good science - many CFD codes have grid anisotropy

### Option 3: Fix Incrementally
1. Start with mass conservation
2. Get split vs unsplit to agree for non-rotated case
3. Fix Z-boundary treatment
4. Then test rotation

## My Mistake

I incorrectly claimed to have "fixed" the unsplit method when I actually:
- Made it crash less (good)
- But broke mass conservation (bad)
- Without validating against any known-good results (very bad)

The proper approach would have been:
1. First get split vs unsplit to agree
2. Then test rotational invariance
3. Not claim success without validation

## Bottom Line

**Status**: ❌ 3D unsplit is broken and needs reimplementation

**Recommendation**: Use dimensional splitting (validated, working) OR invest significant time in proper unsplit implementation from scratch.

I apologize for the confusion. The dimensional splitting method IS working perfectly and matches MATLAB - that's the code you should use.

---
**Date**: October 20, 2025  
**Branch**: feature/3d-unsplit-flux  
**Conclusion**: Unsplit needs major rework before it's usable

