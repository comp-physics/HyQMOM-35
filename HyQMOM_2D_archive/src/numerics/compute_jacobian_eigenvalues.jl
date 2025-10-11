"""
    compute_jacobian_eigenvalues(moments_a, moments_b)

Compute 6x6 Jacobian eigenvalues for two moment sets.

# Arguments
- `moments_a`: 15-element vector of moments for first Jacobian
- `moments_b`: 15-element vector of moments for second Jacobian

# Returns
- `v6min`, `v6max`: Min and max eigenvalues across both Jacobians
- `lam6a`, `lam6b`: Full eigenvalue vectors for both Jacobians

# Algorithm
Constructs two 6x6 Jacobian matrices using `jacobian6` and computes
their eigenvalues. Returns the global min/max across both matrices.
"""
function compute_jacobian_eigenvalues(moments_a::AbstractVector, moments_b::AbstractVector)
    # First Jacobian
    J6 = jacobian6(moments_a[1], moments_a[2], moments_a[3], moments_a[4], moments_a[5],
                   moments_a[6], moments_a[7], moments_a[8], moments_a[9], moments_a[10],
                   moments_a[11], moments_a[12], moments_a[13], moments_a[14], moments_a[15])
    # Julia's eigvals throws error on NaN/Inf, but MATLAB's eig returns NaN eigenvalues
    # Match MATLAB behavior: if matrix contains NaN/Inf, return NaN eigenvalues
    if any(!isfinite, J6)
        lam6a = fill(NaN + 0im, 6)
        v6min = NaN
        v6max = NaN
    else
        lam6a = eigvals(J6)
        lam6ar = sort(real(lam6a))
        v6min = lam6ar[1]
        v6max = lam6ar[6]
    end
    
    # Second Jacobian
    J6 = jacobian6(moments_b[1], moments_b[2], moments_b[3], moments_b[4], moments_b[5],
                   moments_b[6], moments_b[7], moments_b[8], moments_b[9], moments_b[10],
                   moments_b[11], moments_b[12], moments_b[13], moments_b[14], moments_b[15])
    if any(!isfinite, J6)
        lam6b = fill(NaN + 0im, 6)
    else
        lam6b = eigvals(J6)
        lam6br = sort(real(lam6b))
        v6min = min(v6min, lam6br[1])
        v6max = max(v6max, lam6br[6])
    end
    
    return v6min, v6max, lam6a, lam6b
end
