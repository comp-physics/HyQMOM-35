# MPI Domain Decomposition for HyQMOM Solver

This directory contains MPI-style domain decomposition utilities for parallelizing the 2D HyQMOM solver across multiple cores using MATLAB's Parallel Computing Toolbox.

## Overview

The MPI implementation provides:
- **2D Cartesian domain decomposition** splitting the grid into subdomains
- **Efficient halo (ghost cell) exchange** between neighboring processes
- **Automatic load balancing** with block partitioning
- **Flexible boundary conditions** at global domain edges

## Core Components

### 1. Domain Decomposition (`src/setup_mpi_cartesian_2d.m`)

Creates a 2D process grid and assigns subdomains to each worker:

```matlab
decomp = setup_mpi_cartesian_2d(np, halo)
```

**Inputs:**
- `np` - Global grid size (np × np)
- `halo` - Halo width in cells (typically 1)

**Output struct fields:**
- `np_global` - Global grid size
- `halo` - Halo width
- `dims` - Process grid dimensions [Px, Py]
- `coords` - This rank's coordinates [rx, ry] (0-based)
- `rank` - This lab index (1-based)
- `neighbors` - Struct with fields: left, right, down, up (lab indices or -1)
- `local_size` - Interior size [nx_local, ny_local] (without halos)
- `istart_iend` - Global index range for x [i0, i1]
- `jstart_jend` - Global index range for y [j0, j1]

### 2. Halo Exchange (`src/halo_exchange_2d.m`)

Exchanges boundary data between neighboring subdomains:

```matlab
A = halo_exchange_2d(A, decomp, bc)
```

**Inputs:**
- `A` - Local array (nx+2h) × (ny+2h) × nv
- `decomp` - Decomposition struct
- `bc` - Boundary condition struct (optional, default: `struct('type','copy')`)

**Notes:**
- Exchanges left/right first, then up/down
- Fills corners implicitly
- Applies physical BCs at global boundaries

### 3. Physical Boundary Conditions (`src/apply_physical_bc_2d.m`)

Fills halos at global domain boundaries:

```matlab
A = apply_physical_bc_2d(A, decomp, bc)
```

**Supported bc.type:**
- `'copy'` - Neumann-like (copy nearest interior cell)

## Quick Start

### Running the Example

Test the MPI implementation with a simple 5-point stencil:

```matlab
% Start MATLAB, then run:
example_mpi_solver_loop()              % Default: 32×32 grid, 4 workers
example_mpi_solver_loop(64, 10, 3, 4)  % Custom: 64×64, 10 steps, 3 vars, 4 workers
```

Expected output:
```
=== MPI Domain Decomposition Example ===
Grid size: 32 x 32
Workers: 4
Stencil iterations: 5
Variables: 2

Process grid: 2 x 2
Local subdomain size: 16 x 16 (per worker)

Computing serial reference solution...

=== VERIFICATION RESULTS ===
Max absolute error: 1.332e-15
✓ PASS: MPI domain decomposition working correctly!
```

### Integration into Existing Solver

To use MPI domain decomposition in your solver:

1. **Start a parallel pool:**
```matlab
parpool('local', 4);  % 4 workers
```

2. **Wrap your solver in spmd:**
```matlab
spmd
    % Setup decomposition
    decomp = setup_mpi_cartesian_2d(Np, 1);  % halo=1
    nx = decomp.local_size(1);
    ny = decomp.local_size(2);
    
    % Allocate local arrays with halos
    M_local = zeros(nx+2, ny+2, Nmom);
    
    % Scatter global initial conditions to local subdomains
    i0i1 = decomp.istart_iend;
    j0j1 = decomp.jstart_jend;
    M_local(2:nx+1, 2:ny+1, :) = M_global(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :);
    
    % Time stepping loop
    for step = 1:nsteps
        % Exchange halos before flux computation
        M_local = halo_exchange_2d(M_local, decomp);
        
        % Compute fluxes using local data with halos
        % ... your flux computation here ...
        
        % Update interior cells only
        % M_local(2:nx+1, 2:ny+1, :) = ...
    end
    
    % Gather results back to rank 1
    % ... gather logic ...
end
```

3. **Key considerations:**
   - Work on `M_local` instead of global `M`
   - Exchange halos before any operation needing neighbor data
   - Update only interior cells: `M_local(halo+1:halo+nx, halo+1:halo+ny, :)`
   - Handle global boundary conditions automatically via `apply_physical_bc_2d`

## Architecture

### Process Grid Layout

For `numlabs = 4`, the code creates a 2×2 process grid:

```
+-------+-------+
| Rank1 | Rank2 |  ← (coords: [0,1], [1,1])
+-------+-------+
| Rank0 | Rank3 |  ← (coords: [0,0], [1,0])
+-------+-------+
```

### Subdomain Layout with Halos (halo=1)

Each subdomain has interior cells plus ghost cells:

```
+---+-------+---+
| X | Ghost | X |  ← Top halo
+---+-------+---+
| G |       | G |
| h | Inter | h |  ← Left/Right halos
| o | ior   | o |
| s |       | s |
| t |       | t |
+---+-------+---+
| X | Ghost | X |  ← Bottom halo
+---+-------+---+
```

- **Interior:** `A(h+1:h+nx, h+1:h+ny, :)`
- **Left halo:** `A(1:h, h+1:h+ny, :)`
- **Right halo:** `A(h+nx+1:h+nx+h, h+1:h+ny, :)`
- **Bottom halo:** `A(h+1:h+nx, 1:h, :)`
- **Top halo:** `A(h+1:h+nx, h+ny+1:h+ny+h, :)`

### Communication Pattern

Halo exchange occurs in two phases:

1. **Phase 1 (Left/Right):** Exchange vertical strips
2. **Phase 2 (Up/Down):** Exchange horizontal strips (including updated corners)

This ensures corners are filled correctly after both phases.

## Performance Notes

- **No multithreading:** All `parfor` loops have been replaced with `for` loops in `main.m`
- **MPI handles parallelism:** Use `spmd` blocks for parallel execution
- **Load balancing:** Block decomposition distributes remainder cells evenly
- **Communication overhead:** Minimize halo exchanges (only when needed)

## Files Created

```
src/
├── setup_mpi_cartesian_2d.m     # Domain decomposition setup
├── halo_exchange_2d.m            # Halo exchange routine
└── apply_physical_bc_2d.m        # Physical boundary conditions

example_mpi_solver_loop.m         # Example with verification
MPI_README.md                     # This file
```

## Verification

The example script `example_mpi_solver_loop.m` verifies correctness by:
1. Running a 5-point stencil update on decomposed grid with MPI
2. Running the same stencil on the full grid serially
3. Comparing results (should match to machine precision ~10⁻¹⁵)

## Troubleshooting

**Error: "This function must run inside an spmd block"**
- Wrap your MPI code in `spmd ... end`

**Error: No parallel pool**
- Run `parpool('local', N)` before using `spmd`

**High errors in verification**
- Check halo exchange is called before flux computation
- Verify boundary conditions are consistent
- Ensure only interior cells are updated

**Performance not scaling**
- Grid too small for number of workers
- Communication overhead dominates (try larger grid)
- Check load balance with `decomp.local_size`

## Next Steps

For production use:
1. Integrate into full solver (see integration guide above)
2. Add timing/profiling to measure speedup
3. Test with different worker counts
4. Optimize communication patterns if needed
5. Consider adaptive time stepping across ranks

## References

- MATLAB Parallel Computing Toolbox: `spmd`, `labSendReceive`, `labBroadcast`
- HLL Riemann solver: `pas_HLL.m`
- Original solver: `main.m`

