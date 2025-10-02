# MPI Implementation - Current Status

## 🎯 Major Achievement
**✅ Serial vs MPI-1: BITWISE IDENTICAL (0.000e+00)**

This is rock-solid and validates the entire approach!

## 📊 Current State  
**❌ Multi-rank MPI still has errors**

| Configuration | Error vs Serial |
|--------------|-----------------|
| MPI-1 | 0.000e+00 ✅ |
| MPI-2 | 1.140 ❌ |
| MPI-4 | 0.879 ❌ |

## 🔬 What We've Implemented

### Core Infrastructure ✅
1. 2D Cartesian domain decomposition
2. Halo exchange for M, Fx, Fy
3. Wave speed computation in halos
4. Fx/Fy recomputation in halos
5. BC flags in pas_HLL

### Conditional Halo Inclusion ✅ (logic correct, but not working)
```matlab
if has_left_neighbor:
    pass M(halo:halo+nx) to pas_HLL  % [left_halo, interior]
    set apply_bc_left=false
```

The indexing analysis shows this SHOULD work, but results prove it doesn't.

## 🐛 The Mystery

**At boundary cell (12,12):**
- Serial: 0.2051
- MPI-2:  0.1873  
- Difference: **18%!**

This is at the exact processor boundary, suggesting:
1. pas_HLL not using halo data correctly even with BC flags
2. OR: some other issue in how we're calling/extracting from pas_HLL

## 💭 Hypothesis

The BC flags in pas_HLL might not be sufficient. When `apply_bc_left=false`,
we skip applying BCs, but pas_HLL's internal stencil logic might still
have issues at index 1 of the array.

Looking at pas_HLL:
```matlab
for j = 2:Np-1
    lleft(j) = min([vpmin(j), vpmin(j+1)]);
    lright(j) = max([vpmax(j), vpmax(j+1)]);
    Wstar(j,:) = ...
end
```

This loop starts at j=2, so j=1 is NEVER computed for Wstar!
Then later: `Wstar(1,:) = Wstar(2,:)` - this is a BC!

So even with `apply_bc_left=false`, Wstar(1) is still a copy of Wstar(2).

## 🚀 Potential Solution

pas_HLL's loop needs to be `for j = 1:Np-1` when there's valid neighbor
data at index 1. The BC flags should control whether to COMPUTE at the
boundaries, not just whether to apply BCs after.

OR: We need a completely different pas_HLL variant for MPI that handles
processor boundaries explicitly.

## 📁 Files
- `main_mpi.m`: Conditional halo inclusion implemented
- `src/pas_HLL.m`: BC flags added (but insufficient)
- `src/halo_exchange_2d.m`: Working correctly
- Documentation: ANALYSIS.md, MPI_STATUS_FINAL.md, MPI_DEBUGGING_STATUS.md

