# Deep Analysis of MPI Boundary Issue

## The Problem
Serial vs MPI-1: ✅ Perfect (0.000e+00)
MPI ranks differ: ❌ Large error (~0.047 to ~0.796)

## What's happening in Serial

For a 24x24 grid, serial calls `pas_HLL` with arrays of size 24.

At interior cell i=12:
- pas_HLL computes flux between cells 11-12 and 12-13
- Uses M(11), M(12), M(13)
- Uses F(11), F(12), F(13)
- Uses vpmin(11), vpmin(12), vpmin(13)

## What's happening in MPI-2

Grid split: Rank 0 has i=1:12, Rank 1 has i=13:24

**Rank 0, cell 12 (last interior):**
- Local index: nx=12
- pas_HLL gets interior-only array of size 12
- At local index 12:
  - flux_left: between local 11-12 (global 11-12) ✓
  - flux_right: between local 12 and ??? 
  
**The issue:** pas_HLL applies BC at index 12 (the end of the array), 
copying from index 11. But in reality, cell 13 exists on the next rank!

## The Real Solution

pas_HLL MUST see cell 13's data to compute the correct flux at the 12-13 interface.

This means:
1. We can NOT pass interior-only to pas_HLL at processor boundaries
2. We MUST include the halo in the array
3. But we must ensure pas_HLL doesn't apply BCs at processor boundaries

Going back to the BC flags approach, but this time ensuring the data
and indexing are ALL correct.

