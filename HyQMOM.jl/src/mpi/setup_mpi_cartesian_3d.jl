"""
    setup_mpi_cartesian_3d(Np::Int, Nz::Int, halo::Int, comm::MPI.Comm)

Create a 2D Cartesian domain decomposition over MPI ranks (no z decomposition).

# Arguments
- `Np`: Global grid size along x and y (Np x Np)
- `Nz`: Global grid size along z (no decomposition, all ranks have full z)
- `halo`: Halo width in cells (typically 2)
- `comm`: MPI communicator

# Returns
A named tuple with fields:
- `np_global`: Global grid size in x and y
- `nz_global`: Global grid size in z
- `halo`: Halo width
- `dims`: (Px, Py, 1) process grid dimensions (Pz=1, no z decomposition)
- `coords`: (rx, ry, 0) coordinates of this rank in process grid (0-based)
- `rank`: This rank index (0-based)
- `neighbors`: Named tuple with fields left, right, down, up (rank indices or -1)
- `local_size`: (nx_local, ny_local, nz) interior size (without halos)
- `istart_iend`: (i0, i1) global inclusive interior index range for x
- `jstart_jend`: (j0, j1) global inclusive interior index range for y
- `kstart_kend`: (1, nz) global inclusive interior index range for z (always full)

# Notes
- Uses an approximately square process grid (Px, Py, 1) with Px*Py=nprocs
- Uses block decomposition in x and y with remainder cells assigned to lower coords
- No decomposition in z: all ranks have the full z-dimension (Pz=1)

# Example
```julia
decomp = setup_mpi_cartesian_3d(40, 4, 2, MPI.COMM_WORLD)
# For 4 ranks: creates 2×2×1 process grid
# Each rank gets 20×20×4 interior cells
```
"""
function setup_mpi_cartesian_3d(Np::Int, Nz::Int, halo::Int, comm::MPI.Comm)
    # Discover parallel environment
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    # Choose process grid dimensions (Px, Py) - no z decomposition
    Px, Py = choose_process_grid(nprocs)
    Pz = 1  # No decomposition in z
    
    # Rank coordinates (0-based)
    rx = mod(rank, Px)
    ry = div(rank, Px)
    rz = 0  # Always 0 since Pz=1
    
    # Compute local interior sizes using block decomposition with remainders
    nx_local, i0, i1 = block_partition_1d(Np, Px, rx)
    ny_local, j0, j1 = block_partition_1d(Np, Py, ry)
    
    # Z dimension: all ranks have the full z extent
    nz_local = Nz
    k0 = 1
    k1 = Nz
    
    # Neighbor ranks (0-based), -1 if boundary
    # Only x-y neighbors, no z-neighbors since Pz=1
    left = (rx > 0) ? (ry * Px + (rx-1)) : -1
    right = (rx < Px-1) ? (ry * Px + (rx+1)) : -1
    down = (ry > 0) ? ((ry-1) * Px + rx) : -1
    up = (ry < Py-1) ? ((ry+1) * Px + rx) : -1
    
    neighbors = (left=left, right=right, down=down, up=up)
    
    return (
        np_global = Np,
        nz_global = Nz,
        halo = halo,
        dims = (Px, Py, Pz),
        coords = (rx, ry, rz),
        rank = rank,
        neighbors = neighbors,
        local_size = (nx_local, ny_local, nz_local),
        istart_iend = (i0, i1),
        jstart_jend = (j0, j1),
        kstart_kend = (k0, k1)
    )
end

