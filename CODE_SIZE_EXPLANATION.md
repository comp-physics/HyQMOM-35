# Why is main_mpi.m Twice as Long as main.m?

## Line Count
- `main.m`: 341 lines (serial)
- `main_mpi.m`: 645 lines (MPI parallel)
- **Difference: +304 lines (~89% more)**

## Breakdown of Extra Lines

### 1. Halo Wave Speed/Flux Computation (~95 lines)
```matlab
% Lines 210-305: Four nested loops computing in halo cells
for i = 1:halo                      % Left halo
for i = halo+nx+1:nx+2*halo        % Right halo  
for j = 1:halo                      % Bottom halo
for j = halo+ny+1:ny+2*halo        % Top halo
```
**Why needed:** `pas_HLL` stencil requires neighbor data. Must compute wave speeds and fluxes in halo cells after receiving neighbor M data.

### 2. X-Direction Conditional Slicing (~90 lines)
```matlab
% Lines 306-395: Different array slicing based on boundary type
if has_left_neighbor && has_right_neighbor
    % Include halos on both sides
elseif has_left_neighbor
    % Include left halo only
elseif has_right_neighbor  
    % Include right halo only
else
    % Interior only (single rank)
end
```
**Why needed:** `pas_HLL` must know whether to use neighbor data (processor boundary) or apply physical BC (domain boundary).

### 3. Y-Direction Conditional Slicing (~90 lines)
```matlab
% Lines 450-540: Same logic as X-direction for Y-sweeps
```
**Why needed:** Same reason as X-direction, duplicated for Y-sweeps.

### 4. MPI Infrastructure (~50 lines)
- Domain decomposition setup
- Halo exchange calls (M, Fx, Fy)
- Result gathering from all ranks
- Plotting in distributed context

### 5. SPMD Block Overhead (~30 lines)
- Redefining constants inside `spmd` (workers don't inherit workspace)
- Path setup for workers
- Composite array handling

### 6. Additional Code (~50 lines)
- Grid validation for MPI
- Process grid selection
- Parallel pool management

## Why main.m is Simpler

In serial mode:
- ✅ No halos → No halo computation loops
- ✅ No domain boundaries → No conditional slicing
- ✅ No MPI → No communication or gathering
- ✅ Direct array access → No SPMD overhead

Serial `pas_HLL` call (one line):
```matlab
Mnpx(i, :, :) = pas_HLL(squeeze(M(i, :, :)), squeeze(Fx(i, :, :)), ...);
```

MPI `pas_HLL` call (20+ lines per direction):
```matlab
if has_left_neighbor && has_right_neighbor
    i_start = halo; i_end = halo + nx + 1;
    MOM = squeeze(M(i_start:i_end, jh, :));
    FX = squeeze(Fx(i_start:i_end, jh, :));
    vpx_min = vpxmin_ext(vp_start:vp_end, j);
    MNP = pas_HLL(MOM, FX, dt, dx, vpx_min, vpx_max, false, false);
    Mnpx(halo+1:halo+nx, jh, :) = MNP(2:end-1, :);
elseif ...  % 3 more cases
end
```

## Can It Be Simplified?

**No.** Every extra line is necessary for correctness:

1. ❌ Can't remove halo loops → `pas_HLL` stencil needs neighbor wave speeds
2. ❌ Can't remove conditionals → Must distinguish processor vs physical boundaries  
3. ❌ Can't simplify slicing → Array indexing must match boundary type
4. ❌ Can't remove MPI infrastructure → That's the whole point!

## Conclusion

The **2x line count is the inherent cost of distributed parallelism**:
- Halo management
- Boundary type handling  
- Data communication
- Distributed coordination

This is **normal and expected** for MPI implementations. The serial version is simpler because it doesn't deal with domain decomposition.

## Alternative: Remove Serial main.m?

**Recommendation: Keep both**

- `main.m` provides clean serial interface for testing
- Useful for small problems where MPI overhead isn't worth it
- Easier to debug and understand
- Good reference implementation

The unified interface `main(Np, tmax, enable_plots, save_output, use_mpi, num_workers)` lets users choose, so both have value.
