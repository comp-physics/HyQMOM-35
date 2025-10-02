# Halo Exchange Fundamentals

## Core Principle
A halo (ghost cell) region contains a COPY of data from neighboring ranks.
After halo exchange, the overlapping regions must be BITWISE IDENTICAL.

## What Must Be Exchanged

For a finite volume method with flux computation:

1. **Conservative variables (M)**: PRIMARY data
   - Each rank owns interior cells
   - Halos contain copies of neighbors' interior cells
   - MUST be exchanged before any operation that uses neighbors

2. **Fluxes (Fx, Fy)**: DERIVED from M
   - Computed locally from M at each cell
   - If flux at cell boundary depends on M from both sides, halos needed
   - Question: Does Fx/Fy computation use neighbor M values?

3. **Wave speeds (vpmin, vpmax)**: DERIVED from M
   - Computed locally from M at each cell
   - Used in HLL flux computation
   - Question: Does wave speed computation use neighbor M values?

## Key Questions for This Code

1. **Does Fx/Fy computation (lines 186-193 of main_mpi.m) use M from neighbor cells?**
   - Looking at code: Fx(ih,jh,:) = Mx from closure of M(ih,jh,:)
   - Answer: NO - it's a pointwise operation on single cell
   - Conclusion: Fx/Fy halos should be computed from exchanged M halos, OR exchanged directly

2. **Does pas_HLL need neighbor Fx/Fy values?**
   - Looking at pas_HLL line 11: uses F(j,:) and F(j+1,:)
   - Answer: YES - it needs F at j and j+1 (neighbor in sweep direction)
   - Conclusion: F halos MUST contain correct neighbor F values

3. **Does pas_HLL need neighbor wave speeds?**
   - Looking at pas_HLL line 6: uses vpmin(j) and vpmin(j+1)
   - Answer: YES - needs wave speeds at j and j+1
   - Conclusion: Wave speed halos MUST contain correct neighbor values

## Current Implementation Check

Need to verify:
- [ ] M halos exchanged before Fx/Fy computation?
- [ ] Fx/Fy halos filled (either computed from M halos OR exchanged)?
- [ ] Wave speed halos filled before pas_HLL?
- [ ] All exchanges happen in correct order?

## Debugging Strategy

1. Increase halo to 2 to make bugs more obvious
2. Trace exactly what pas_HLL sees at processor boundaries
3. Verify F(1) contains neighbor's F value, not BC
4. Check if wave speeds at halos are correct
