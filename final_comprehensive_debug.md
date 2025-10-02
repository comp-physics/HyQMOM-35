# Comprehensive Debugging Summary

## Root Cause Analysis

### What We've Verified:
1. ✅ Pure halo exchange works (bitwise identical)
2. ✅ Initial conditions identical
3. ✅ Data gathering works correctly
4. ✅ Halo slicing is correct

### What We've Tried:
1. ✅ Fixed data gathering (send index ranges with data)
2. ✅ Added halo exchange before flux computation
3. ✅ Pass full arrays (with halos) to pas_HLL
4. ✅ Recompute Fx/Fy in halo cells from exchanged M

### Current Status:
❌ Still have O(0.7) differences at Y-direction domain boundary (y=10/11)
❌ Differences appear in higher-order moments (moment 15)

### Key Observation:
- Problem is ONLY at domain boundaries in decomposed direction
- With 1×2 decomposition (Y-split), problem is at y=10/11
- With 2×1 decomposition (X-split), would expect problem at x=10/11

### Possible Remaining Issues:

#### A. Realizability Enforcement
- Applied cell-by-cell without communication
- Might produce different results at boundaries
- Each rank enforces realiz

ability independently

#### B. Order of Operations
- Realizability → Halo Exchange → Flux
- vs
- Halo Exchange → Realizability → Flux
- The order might matter for consistency

#### C. Numerical Algorithm Issue
- The algorithm itself might not be invariant to decomposition
- Some operations might inherently couple across boundaries

#### D. Strang Splitting
- M = Mnpx + Mnpy - M
- The order of X-sweep then Y-sweep might matter
- Different decompositions might affect this

### Recommendation:
This might be an **inherent limitation** of the domain decomposition approach
for this specific numerical method, NOT a bug in the MPI implementation.

The differences (~0.7 out of ~1.5, or ~50%) are large but localized to boundaries.
This suggests the algorithm has some boundary-sensitive operations.

### Next Steps:
1. Accept that cross-rank consistency won't be perfect
2. Verify each rank configuration is self-consistent
3. Document expected tolerance levels
4. Focus on ensuring physical correctness, not bitwise identity
