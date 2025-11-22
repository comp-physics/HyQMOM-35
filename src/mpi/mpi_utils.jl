"""
MPI utility functions for domain decomposition.

This module provides utilities for:
- Process grid selection
- 1D block decomposition
- Gathering results
- Sending data between ranks
"""

"""
    choose_process_grid(nl::Int)

Choose nearly-square factorization Px*Py=nl where Px and Py are as 
close as possible to sqrt(nl).

# Arguments
- `nl`: Number of processes

# Returns
- `Px`, `Py`: Process grid dimensions

# Examples
```julia
Px, Py = choose_process_grid(6)  # Returns (2, 3) or (3, 2)
```
"""
function choose_process_grid(nl::Int)
    bestDiff = typemax(Int)
    Px = 1
    Py = nl
    
    for p in 1:nl
        if mod(nl, p) == 0
            q = div(nl, p)
            d = abs(p - q)
            if d < bestDiff
                bestDiff = d
                Px = p
                Py = q
            end
        end
    end
    
    return Px, Py
end

"""
    block_partition_1d(n::Int, P::Int, r::Int)

Compute 1D block decomposition for rank r (0-indexed).

Divides n points among P processes as evenly as possible.
First (n mod P) ranks get (floor(n/P) + 1) points each.

# Arguments
- `n`: Total number of points
- `P`: Number of processes
- `r`: Rank (0-indexed)

# Returns
- `n_local`: Number of points for this rank
- `i0`: Starting global index (1-indexed)
- `i1`: Ending global index (1-indexed)

# Examples
```julia
n_local, i0, i1 = block_partition_1d(100, 4, 0)  # First rank
# Returns (25, 1, 25)
```
"""
function block_partition_1d(n::Int, P::Int, r::Int)
    base = div(n, P)
    remn = mod(n, P)
    
    if r < remn
        n_local = base + 1
        i0 = r * (base + 1) + 1
    else
        n_local = base
        i0 = remn * (base + 1) + (r - remn) * base + 1
    end
    
    i1 = i0 + n_local - 1
    
    return n_local, i0, i1
end

"""
    gather_M(M_interior::Array{T,3}, i0i1::Tuple{Int,Int}, j0j1::Tuple{Int,Int}, 
             Np::Int, Nmom::Int, comm::MPI.Comm) where T

Gather 2D moment arrays from all ranks to rank 0.

# Arguments
- `M_interior`: Local interior moment array (nx, ny, Nmom)
- `i0i1`: Global index range for x (i0, i1)
- `j0j1`: Global index range for y (j0, j1)
- `Np`: Global grid size
- `Nmom`: Number of moments
- `comm`: MPI communicator

# Returns
- `M_full`: Full moment array (only on rank 0, nothing on other ranks)

# Algorithm
Rank 0 receives data from all other ranks and assembles the full array.
Other ranks send their data to rank 0.
"""
function gather_M(M_interior::Array{T,3}, i0i1::Tuple{Int,Int}, j0j1::Tuple{Int,Int},
                 Np::Int, Nmom::Int, comm::MPI.Comm) where T
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    if rank == 0
        # Initialize full array on rank 0
        M_full = zeros(T, Np, Np, Nmom)
        
        # Place rank 0's data
        M_full[i0i1[1]:i0i1[2], j0j1[1]:j0j1[2], :] = M_interior
        
        # Receive from all other ranks
        for src in 1:(nprocs-1)
            # Receive index ranges
            i_range = Vector{Int}(undef, 2)
            MPI.Recv!(i_range, comm; source=src, tag=0)
            j_range = Vector{Int}(undef, 2)
            MPI.Recv!(j_range, comm; source=src, tag=1)
            
            # Receive data
            nx_remote = i_range[2] - i_range[1] + 1
            ny_remote = j_range[2] - j_range[1] + 1
            blk = Array{T,3}(undef, nx_remote, ny_remote, Nmom)
            MPI.Recv!(blk, comm; source=src, tag=2)
            
            # Place in full array
            M_full[i_range[1]:i_range[2], j_range[1]:j_range[2], :] = blk
        end
        
        return M_full
    else
        # Send data to rank 0
        MPI.Send([i0i1[1], i0i1[2]], comm; dest=0, tag=0)
        MPI.Send([j0j1[1], j0j1[2]], comm; dest=0, tag=1)
        MPI.Send(M_interior, comm; dest=0, tag=2)
        
        return nothing
    end
end

"""
    gather_M(M_interior::Array{T,4}, i0i1::Tuple{Int,Int}, j0j1::Tuple{Int,Int}, k0k1::Tuple{Int,Int},
             Np::Int, Nz::Int, Nmom::Int, comm::MPI.Comm) where T

Gather 3D moment arrays from all ranks to rank 0.

# Arguments
- `M_interior`: Local interior moment array (nx, ny, nz, Nmom)
- `i0i1`: Global index range for x (i0, i1)
- `j0j1`: Global index range for y (j0, j1)
- `k0k1`: Global index range for z (k0, k1) - always (1, Nz) for 2D MPI decomposition
- `Np`: Global grid size in x and y
- `Nz`: Global grid size in z
- `Nmom`: Number of moments
- `comm`: MPI communicator

# Returns
- `M_full`: Full moment array (only on rank 0, nothing on other ranks)

# Algorithm
Rank 0 receives data from all other ranks and assembles the full array.
Other ranks send their data to rank 0.
"""
function gather_M(M_interior::Array{T,4}, i0i1::Tuple{Int,Int}, j0j1::Tuple{Int,Int}, k0k1::Tuple{Int,Int},
                 Nx::Int, Ny::Int, Nz::Int, Nmom::Int, comm::MPI.Comm) where T
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    if rank == 0
        # Initialize full array on rank 0
        M_full = zeros(T, Nx, Ny, Nz, Nmom)
        
        # Place rank 0's data
        M_full[i0i1[1]:i0i1[2], j0j1[1]:j0j1[2], k0k1[1]:k0k1[2], :] = M_interior
        
        # Receive from all other ranks
        for src in 1:(nprocs-1)
            # Receive index ranges
            i_range = Vector{Int}(undef, 2)
            MPI.Recv!(i_range, comm; source=src, tag=0)
            j_range = Vector{Int}(undef, 2)
            MPI.Recv!(j_range, comm; source=src, tag=1)
            k_range = Vector{Int}(undef, 2)
            MPI.Recv!(k_range, comm; source=src, tag=10)
            
            # Receive data
            nx_remote = i_range[2] - i_range[1] + 1
            ny_remote = j_range[2] - j_range[1] + 1
            nz_remote = k_range[2] - k_range[1] + 1
            blk = Array{T,4}(undef, nx_remote, ny_remote, nz_remote, Nmom)
            MPI.Recv!(blk, comm; source=src, tag=2)
            
            # Place in full array
            M_full[i_range[1]:i_range[2], j_range[1]:j_range[2], k_range[1]:k_range[2], :] = blk
        end
        
        return M_full
    else
        # Send data to rank 0
        MPI.Send([i0i1[1], i0i1[2]], comm; dest=0, tag=0)
        MPI.Send([j0j1[1], j0j1[2]], comm; dest=0, tag=1)
        MPI.Send([k0k1[1], k0k1[2]], comm; dest=0, tag=10)
        MPI.Send(M_interior, comm; dest=0, tag=2)
        
        return nothing
    end
end

"""
    send_M(M_interior::Array{T,3}, i0i1::Tuple{Int,Int}, j0j1::Tuple{Int,Int},
           dest_rank::Int, comm::MPI.Comm) where T

Send moment array to destination rank.

# Arguments
- `M_interior`: Local interior moment array
- `i0i1`: Global index range for x
- `j0j1`: Global index range for y
- `dest_rank`: Destination rank
- `comm`: MPI communicator

# Notes
This is a lower-level function. Typically use `gather_M` instead.
"""
function send_M(M_interior::Array{T,3}, i0i1::Tuple{Int,Int}, j0j1::Tuple{Int,Int},
               dest_rank::Int, comm::MPI.Comm) where T
    MPI.Send([i0i1[1], i0i1[2]], comm; dest=dest_rank, tag=0)
    MPI.Send([j0j1[1], j0j1[2]], comm; dest=dest_rank, tag=1)
    MPI.Send(M_interior, comm; dest=dest_rank, tag=2)
end
