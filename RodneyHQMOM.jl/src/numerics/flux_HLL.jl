"""
    flux_HLL(Wstar, W, l1, l2, F, N)

Compute HLL flux from intermediate states and wave speeds.

This is a simple 5-line function that computes the numerical flux
using the HLL (Harten-Lax-van Leer) scheme.

# Arguments
- `Wstar`: Intermediate state array
- `W`: State array
- `l1`: Left wave speeds
- `l2`: Right wave speeds
- `F`: Flux array
- `N`: Number of cells

# Returns
- `Flux`: Numerical flux array
"""
function flux_HLL(Wstar, W, l1, l2, F, N)
    Flux = zeros(size(F, 1)-1, size(F, 2))
    
    for j = 1:N-1
        for k = 1:size(F, 2)
            Flux[j,k] = 0.5*(F[j,k] + F[j+1,k]) - 
                        0.5*(abs(l1[j])*(Wstar[j,k] - W[j,k]) - 
                             abs(l2[j])*(Wstar[j,k] - W[j+1,k]))
        end
    end
    
    return Flux
end
