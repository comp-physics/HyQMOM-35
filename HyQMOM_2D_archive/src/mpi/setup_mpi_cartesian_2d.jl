"""
    setup_mpi_cartesian_2d(np::Int, halo::Int, comm::MPI.Comm)

Create a 2D Cartesian domain decomposition.

# Arguments
- `np`: Global grid size along x and y (np x np)
- `halo`: Halo width in cells (typically 1)
- `comm`: MPI communicator

# Returns
A named tuple with fields:
- `np_global`: Scalar, global grid size
- `halo`: Scalar, halo width
- `dims`: (Px, Py) process grid dimensions
- `coords`: (rx, ry) coordinates of this rank in process grid (0-based)
- `rank`: This rank index (0-based)
- `size`: Total number of ranks
- `neighbors`: Named tuple with fields: left, right, down, up (rank indices or -1)
- `local_size`: (nx_local, ny_local) interior size (without halos)
- `istart_iend`: (i0, i1) global inclusive interior index range for x
- `jstart_jend`: (j0, j1) global inclusive interior index range for y

# Algorithm
1. Choose approximately square process grid [Px, Py] with Px*Py = nprocs
2. Use block decomposition with remainder cells assigned to lower coordinates
3. Compute neighbor ranks using Cartesian topology

# Examples
```julia
using MPI
MPI.Init()
comm = MPI.COMM_WORLD
decomp = setup_mpi_cartesian_2d(120, 1, comm)
println("Rank \$(decomp.rank) has local size \$(decomp.local_size)")
```

# Notes
- Uses 0-based rank indexing (Julia convention for MPI)
- Neighbor ranks are -1 at physical boundaries
- Process grid is chosen to be as square as possible
"""
function setup_mpi_cartesian_2d(np::Int, halo::Int, comm::MPI.Comm)
    # Discover parallel environment
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    # Choose process grid dimensions (Px, Py)
    Px, Py = choose_process_grid(nprocs)
    
    # Rank coordinates (0-based)
    rx = mod(rank, Px)
    ry = div(rank, Px)
    
    # Compute local interior sizes using block decomposition with remainders
    nx_local, i0, i1 = block_partition_1d(np, Px, rx)
    ny_local, j0, j1 = block_partition_1d(np, Py, ry)
    
    # Neighbor ranks (0-based), -1 if boundary
    left_neighbor = (rx > 0) ? rank_from_coords(rx-1, ry, Px) : -1
    right_neighbor = (rx < Px-1) ? rank_from_coords(rx+1, ry, Px) : -1
    down_neighbor = (ry > 0) ? rank_from_coords(rx, ry-1, Px) : -1
    up_neighbor = (ry < Py-1) ? rank_from_coords(rx, ry+1, Px) : -1
    
    neighbors = (left=left_neighbor, right=right_neighbor, 
                 down=down_neighbor, up=up_neighbor)
    
    return (np_global = np,
            halo = halo,
            dims = (Px, Py),
            coords = (rx, ry),
            rank = rank,
            size = nprocs,
            neighbors = neighbors,
            local_size = (nx_local, ny_local),
            istart_iend = (i0, i1),
            jstart_jend = (j0, j1))
end

"""
    rank_from_coords(rx::Int, ry::Int, Px::Int)

Convert (rx, ry) 0-based coordinates to 0-based rank index.

# Arguments
- `rx`: X-coordinate in process grid (0-based)
- `ry`: Y-coordinate in process grid (0-based)
- `Px`: Number of processes in X-direction

# Returns
- Rank index (0-based)

# Algorithm
Uses row-major ordering: rank = ry * Px + rx
"""
function rank_from_coords(rx::Int, ry::Int, Px::Int)
    return ry * Px + rx
end
