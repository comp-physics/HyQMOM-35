"""
    halo_exchange_2d!(A::Array{T,3}, decomp, bc=nothing) where T

Exchange halos for a 2D subdomain with nv variables.

# Arguments
- `A`: Local array (nx+2h) × (ny+2h) × nv, interior: A[h+1:h+nx, h+1:h+ny, :]
- `decomp`: Struct from `setup_mpi_cartesian_2d`
- `bc`: (optional) Boundary condition type, default: `:copy`

# Returns
- Updates `A` in-place with exchanged halos

# Algorithm
1. Apply physical BC at global boundaries before exchange
2. Perform left/right exchange first
3. Then perform up/down exchange
4. Corners are filled implicitly after second phase

# Notes
- Uses blocking send/receive (MPI.Send! and MPI.Recv!)
- For single-rank case, only applies physical BCs
- Physical boundaries use copy BC by default

# Examples
```julia
A = zeros(nx+2*halo, ny+2*halo, Nmom)
# ... fill interior ...
halo_exchange_2d!(A, decomp)
# Now A has valid halo data
```
"""
function halo_exchange_2d!(A::Array{T,3}, decomp, bc=:copy) where T
    h = decomp.halo
    nx = decomp.local_size[1]
    ny = decomp.local_size[2]
    
    if h == 0
        return A
    end
    
    # Apply physical BC to halos at global boundaries before exchange
    apply_physical_bc_2d!(A, decomp, bc)
    
    if decomp.size == 1
        return A
    end
    
    comm = MPI.COMM_WORLD
    
    # X-direction exchange
    left_neighbor = decomp.neighbors.left
    right_neighbor = decomp.neighbors.right
    
    # Send interior boundary data to neighbors
    if left_neighbor != -1
        send_buf = A[h+1:h+h, h+1:h+ny, :]
        MPI.Send(send_buf, left_neighbor, 0, comm)
    end
    
    if right_neighbor != -1
        send_buf = A[h+nx-h+1:h+nx, h+1:h+ny, :]
        MPI.Send(send_buf, right_neighbor, 1, comm)
    end
    
    # Receive halo data from neighbors
    if left_neighbor != -1
        recv_buf = similar(A, h, ny, size(A, 3))
        MPI.Recv!(recv_buf, left_neighbor, 1, comm)
        A[1:h, h+1:h+ny, :] = recv_buf
    end
    
    if right_neighbor != -1
        recv_buf = similar(A, h, ny, size(A, 3))
        MPI.Recv!(recv_buf, right_neighbor, 0, comm)
        A[h+nx+1:h+nx+h, h+1:h+ny, :] = recv_buf
    end
    
    # Y-direction exchange
    down_neighbor = decomp.neighbors.down
    up_neighbor = decomp.neighbors.up
    
    # Send interior boundary data to neighbors
    if down_neighbor != -1
        send_buf = A[h+1:h+nx, h+1:h+h, :]
        MPI.Send(send_buf, down_neighbor, 2, comm)
    end
    
    if up_neighbor != -1
        send_buf = A[h+1:h+nx, h+ny-h+1:h+ny, :]
        MPI.Send(send_buf, up_neighbor, 3, comm)
    end
    
    # Receive halo data from neighbors
    if down_neighbor != -1
        recv_buf = similar(A, nx, h, size(A, 3))
        MPI.Recv!(recv_buf, down_neighbor, 3, comm)
        A[h+1:h+nx, 1:h, :] = recv_buf
    end
    
    if up_neighbor != -1
        recv_buf = similar(A, nx, h, size(A, 3))
        MPI.Recv!(recv_buf, up_neighbor, 2, comm)
        A[h+1:h+nx, h+ny+1:h+ny+h, :] = recv_buf
    end
    
    return A
end

"""
    apply_physical_bc_2d!(A::Array{T,3}, decomp, bc) where T

Apply physical boundary conditions at global boundaries.

# Arguments
- `A`: Local array with halos
- `decomp`: Decomposition structure
- `bc`: Boundary condition type (`:copy`, `:periodic`, etc.)

# Algorithm
Only applies BCs at physical (global) boundaries, not at processor boundaries.
Uses copy BC: halo cells copy from nearest interior cell.

# Notes
This is called before halo exchange to set physical boundary values.
"""
function apply_physical_bc_2d!(A::Array{T,3}, decomp, bc) where T
    h = decomp.halo
    nx = decomp.local_size[1]
    ny = decomp.local_size[2]
    
    # Copy BC: halo cells copy from nearest interior cell
    if bc == :copy
        # Left boundary (global)
        if decomp.neighbors.left == -1
            for i in 1:h
                A[i, h+1:h+ny, :] = A[h+1:h+1, h+1:h+ny, :]
            end
        end
        
        # Right boundary (global)
        if decomp.neighbors.right == -1
            for i in 1:h
                A[h+nx+i, h+1:h+ny, :] = A[h+nx:h+nx, h+1:h+ny, :]
            end
        end
        
        # Down boundary (global)
        if decomp.neighbors.down == -1
            for j in 1:h
                A[h+1:h+nx, j, :] = A[h+1:h+nx, h+1:h+1, :]
            end
        end
        
        # Up boundary (global)
        if decomp.neighbors.up == -1
            for j in 1:h
                A[h+1:h+nx, h+ny+j, :] = A[h+1:h+nx, h+ny:h+ny, :]
            end
        end
        
        # Corners (copy from nearest interior corner)
        if decomp.neighbors.left == -1 && decomp.neighbors.down == -1
            A[1:h, 1:h, :] .= A[h+1:h+1, h+1:h+1, :]
        end
        if decomp.neighbors.right == -1 && decomp.neighbors.down == -1
            A[h+nx+1:h+nx+h, 1:h, :] .= A[h+nx:h+nx, h+1:h+1, :]
        end
        if decomp.neighbors.left == -1 && decomp.neighbors.up == -1
            A[1:h, h+ny+1:h+ny+h, :] .= A[h+1:h+1, h+ny:h+ny, :]
        end
        if decomp.neighbors.right == -1 && decomp.neighbors.up == -1
            A[h+nx+1:h+nx+h, h+ny+1:h+ny+h, :] .= A[h+nx:h+nx, h+ny:h+ny, :]
        end
    end
    
    return A
end
