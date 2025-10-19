"""
Full 35-moment tensor rotation utilities for validation tests.

These functions transform ALL 35 moments under coordinate rotations,
including proper tensor transformations for stress, flux, and higher moments.
"""

"""
    rotate_moments_z90(M::AbstractVector) -> Vector

Apply Z-axis rotation (+90°) to full 35-moment vector.

Coordinate transformation: (x, y, z) → (y, -x, z)
Velocity transformation:   (u, v, w) → (v, -u, w)

All moments M_{ijk} = ⟨u^i v^j w^k⟩ transform via:
    M_{ijk}^new = ⟨(v)^i (-u)^j (w)^k⟩
                = (-1)^j ⟨v^i u^j w^k⟩
                = (-1)^j M_{jik}^old

# Arguments
- `M`: 35-element moment vector in standard order

# Returns
- `M_rot`: 35-element rotated moment vector

# Moment ordering (i,j,k):
```
 1: (0,0,0)   6: (0,1,0)  11: (0,2,1)  16: (0,0,1)  21: (0,1,2)  26: (0,1,1)  31: (0,2,2)
 2: (1,0,0)   7: (1,1,0)  12: (0,3,1)  17: (1,0,1)  22: (1,2,1)  27: (1,0,2)  32: (0,1,2)
 3: (2,0,0)   8: (2,1,0)  13: (0,3,0)  18: (2,0,1)  23: (0,0,3)  28: (2,0,2)  33: (1,1,2)
 4: (3,0,0)   9: (3,1,0)  14: (1,2,0)  19: (3,0,1)  24: (1,0,3)  29: (0,2,2)  34: (0,1,3)
 5: (4,0,0)  10: (0,2,0)  15: (0,4,0)  20: (0,0,2)  25: (0,0,4)  30: (1,1,1)  35: (0,3,2)
```
"""
function rotate_moments_z90(M::AbstractVector)
    M_rot = similar(M)
    
    # Systematically apply: M_ijk^new = (-1)^j M_jik^old
    # Using ordering from InitializeM4_35.jl
    
    M_rot[1] = M[1]    # M000 → M000
    M_rot[2] = -M[6]   # M100 → -M010
    M_rot[3] = M[10]   # M200 → M020
    M_rot[4] = -M[13]  # M300 → -M030
    M_rot[5] = M[15]   # M400 → M040
    M_rot[6] = -M[2]   # M010 → -M100
    M_rot[7] = M[7]    # M110 → M110 (actually  -(-M110) = M110)
    M_rot[8] = -M[11]  # M210 → -M120
    M_rot[9] = M[14]   # M310 → M130
    M_rot[10] = M[3]   # M020 → M200
    M_rot[11] = -M[8]  # M120 → -M210
    M_rot[12] = M[12]  # M220 → M220
    M_rot[13] = -M[4]  # M030 → -M300
    M_rot[14] = M[9]   # M130 → M310
    M_rot[15] = M[5]   # M040 → M400
    M_rot[16] = M[16]  # M001 → M001
    M_rot[17] = -M[26] # M101 → -M011
    M_rot[18] = M[29]  # M201 → M021
    M_rot[19] = -M[31] # M301 → -M031
    M_rot[20] = M[20]  # M002 → M002
    M_rot[21] = -M[32] # M102 → -M012
    M_rot[22] = M[35]  # M202 → M022
    M_rot[23] = M[23]  # M003 → M003
    M_rot[24] = -M[34] # M103 → -M013
    M_rot[25] = M[25]  # M004 → M004
    M_rot[26] = -M[17] # M011 → -M101
    M_rot[27] = M[27]  # M111 → M111
    M_rot[28] = -M[30] # M211 → -M121
    M_rot[29] = M[18]  # M021 → M201
    M_rot[30] = -M[28] # M121 → -M211
    M_rot[31] = -M[19] # M031 → -M301
    M_rot[32] = -M[21] # M012 → -M102
    M_rot[33] = M[33]  # M112 → M112
    M_rot[34] = -M[24] # M013 → -M103
    M_rot[35] = M[22]  # M022 → M202
    
    return M_rot
end

"""
    rotate_moments_z90_inverse(M::AbstractVector) -> Vector

Apply inverse Z-axis rotation (-90°) to full 35-moment vector.

Coordinate transformation: (x, y, z) → (-y, x, z)
Velocity transformation:   (u, v, w) → (-v, u, w)

This is the inverse of `rotate_moments_z90`.
"""
function rotate_moments_z90_inverse(M::AbstractVector)
    M_rot = similar(M)
    
    # Inverse rotation: (u,v,w) → (-v,u,w)
    # So M_ijk^new = (-1)^i M_jik^old (note: i and j swap, sign from i)
    
    M_rot[1] = M[1]    # M000 → M000
    M_rot[2] = -M[6]   # M100 → -M010
    M_rot[3] = M[10]   # M200 → M020
    M_rot[4] = -M[13]  # M300 → -M030
    M_rot[5] = M[15]   # M400 → M040
    M_rot[6] = -M[2]   # M010 → -M100
    M_rot[7] = M[7]    # M110 → M110
    M_rot[8] = -M[11]  # M210 → -M120
    M_rot[9] = M[14]   # M310 → M130
    M_rot[10] = M[3]   # M020 → M200
    M_rot[11] = -M[8]  # M120 → -M210
    M_rot[12] = M[12]  # M220 → M220
    M_rot[13] = -M[4]  # M030 → -M300
    M_rot[14] = M[9]   # M130 → M310
    M_rot[15] = M[5]   # M040 → M400
    M_rot[16] = M[16]  # M001 → M001
    M_rot[17] = -M[26] # M101 → -M011
    M_rot[18] = M[29]  # M201 → M021
    M_rot[19] = -M[31] # M301 → -M031
    M_rot[20] = M[20]  # M002 → M002
    M_rot[21] = -M[32] # M102 → -M012
    M_rot[22] = M[35]  # M202 → M022
    M_rot[23] = M[23]  # M003 → M003
    M_rot[24] = -M[34] # M103 → -M013
    M_rot[25] = M[25]  # M004 → M004
    M_rot[26] = -M[17] # M011 → -M101
    M_rot[27] = M[27]  # M111 → M111
    M_rot[28] = -M[30] # M211 → -M121
    M_rot[29] = M[18]  # M021 → M201
    M_rot[30] = -M[28] # M121 → -M211
    M_rot[31] = -M[19] # M031 → -M301
    M_rot[32] = -M[21] # M012 → -M102
    M_rot[33] = M[33]  # M112 → M112
    M_rot[34] = -M[24] # M013 → -M103
    M_rot[35] = M[22]  # M022 → M202
    
    return M_rot
end

"""
    permute_grid_z90(M::Array{T,4}, Nx::Int, Ny::Int, Nz::Int) -> Array{T,4}

Apply Z-rotation (+90°) grid permutation to 4D moment array.

Grid transformation: (i, j, k) → (j, Nx-i+1, k)

This handles the spatial rearrangement of grid cells under rotation.
After this permutation, call `rotate_moments_z90` on each cell's moment vector.

# Arguments
- `M`: 4D array (Nx, Ny, Nz, Nmom)
- `Nx, Ny, Nz`: Grid dimensions

# Returns
- `M_perm`: Permuted 4D array matching rotated coordinate system
"""
function permute_grid_z90(M::Array{T,4}, Nx::Int, Ny::Int, Nz::Int) where T
    Nmom = size(M, 4)
    M_perm = zeros(T, Nx, Ny, Nz, Nmom)
    
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        i_new = j
        j_new = Nx - i + 1
        k_new = k
        
        M_perm[i_new, j_new, k_new, :] = M[i, j, k, :]
    end
    
    return M_perm
end

"""
    permute_grid_z90_inverse(M::Array{T,4}, Nx::Int, Ny::Int, Nz::Int) -> Array{T,4}

Apply inverse Z-rotation (-90°) grid permutation.

Grid transformation: (i, j, k) → (Ny-j+1, i, k)
"""
function permute_grid_z90_inverse(M::Array{T,4}, Nx::Int, Ny::Int, Nz::Int) where T
    Nmom = size(M, 4)
    M_perm = zeros(T, Nx, Ny, Nz, Nmom)
    
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        i_new = Ny - j + 1
        j_new = i
        k_new = k
        
        M_perm[i_new, j_new, k_new, :] = M[i, j, k, :]
    end
    
    return M_perm
end

"""
    rotate_full_state_z90_inverse(M::Array{T,4}) -> Array{T,4}

Complete state rotation: grid permutation + moment transformation.

This applies BOTH:
1. Grid permutation (spatial rearrangement)
2. Moment tensor rotation (velocity/tensor transformation)

Use this for comparing rotated simulation results back to reference frame.
"""
function rotate_full_state_z90_inverse(M::Array{T,4}) where T
    Nx, Ny, Nz, Nmom = size(M)
    
    # Step 1: Permute grid
    M_perm = permute_grid_z90_inverse(M, Nx, Ny, Nz)
    
    # Step 2: Rotate moments at each grid point
    M_rot = similar(M_perm)
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        M_rot[i, j, k, :] = rotate_moments_z90_inverse(M_perm[i, j, k, :])
    end
    
    return M_rot
end

