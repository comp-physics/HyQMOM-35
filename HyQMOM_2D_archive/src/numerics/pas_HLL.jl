"""
    pas_HLL(M, F, dt, dx, vpmin, vpmax; apply_bc_left=true, apply_bc_right=true)

HLL flux update scheme.

Computes the updated moments using the HLL (Harten-Lax-van Leer) flux scheme.

# Arguments
- `M`: Moment array (Np x Nmom)
- `F`: Flux array (Np x Nmom)
- `dt`: Time step
- `dx`: Spatial step
- `vpmin`: Minimum eigenvalues (Np)
- `vpmax`: Maximum eigenvalues (Np)
- `apply_bc_left`: Apply boundary conditions at left boundary (default: true)
- `apply_bc_right`: Apply boundary conditions at right boundary (default: true)

# Returns
- `Mp`: Updated moment array

# Notes
- Set `apply_bc_left=false` or `apply_bc_right=false` for processor boundaries in MPI
"""
function pas_HLL(M, F, dt, dx, vpmin, vpmax; apply_bc_left=true, apply_bc_right=true)
    Np = size(M, 1)
    Nmom = size(M, 2)
    Wstar = zeros(Np, Nmom)
    lleft = zeros(Np)
    lright = zeros(Np)
    
    # Determine loop range based on boundary conditions
    # If left processor boundary (apply_bc_left=false), extend loop to START at 1
    j_start = 2
    if !apply_bc_left
        j_start = 1  # Start at 1 to compute Wstar(1) using left neighbor M(1)
    end
    
    # For right boundary, we always end at Np-1 because computing at Np would need M(Np+1)
    j_end = Np - 1
    
    # Compute Wstar using stencil (needs M(j) and M(j+1))
    for j = j_start:j_end
        lleft[j] = min(vpmin[j], vpmin[j+1])
        lright[j] = max(vpmax[j], vpmax[j+1])
        
        # Wstar
        if abs(lleft[j] - lright[j]) > 1e-10
            Wstar[j,:] = (lleft[j]*M[j,:] - lright[j]*M[j+1,:]) / (lleft[j] - lright[j]) -
                         (F[j,:] - F[j+1,:]) / (lleft[j] - lright[j])
        else
            Wstar[j,:] .= 0.0
        end
    end
    
    # Apply boundary conditions ONLY at physical boundaries
    if apply_bc_left
        Wstar[1,:] = Wstar[2,:]
    end
    if apply_bc_right
        Wstar[Np,:] = Wstar[Np-1,:]
    end
    
    if apply_bc_left
        F[1,:] = F[2,:]
    end
    if apply_bc_right
        F[Np,:] = F[Np-1,:]
    end
    
    Flux = flux_HLL(Wstar, M, lleft, lright, F, Np)
    
    Mp = copy(M)
    Mp[2:Np-1,:] = M[2:Np-1,:] - dt/dx * (Flux[2:Np-1,:] - Flux[1:Np-2,:])
    
    # Apply solution BCs only at physical boundaries
    if apply_bc_left
        Mp[1,:] = Mp[2,:]
    end
    if apply_bc_right
        Mp[Np,:] = Mp[Np-1,:]
    end
    
    return Mp
end
