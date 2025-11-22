"""
    realizability_S210(S110, S101, S011, S300, S210, S201, H200, beta)

Check and correct realizability of S210 and S201 moments.

# Arguments
- S110, S101, S011: Second-order standardized moments
- S300: Third-order moment along first axis
- S210, S201: Third-order cross-moments to correct
- H200: Variance-related quantity (H200 = max(eps, S400 - S300^2 - 1))
- beta: Scaling factor (typically 1.0)

# Returns
- S210r, S201r: Corrected moments

# Algorithm
Uses matrix square root and quadratic form to determine feasible region.
Scales the deviation from the mean by xr in [0,1] to satisfy realizability.
"""
@fastmath function realizability_S210(S110, S101, S011, S300, S210, S201, H200, beta)
    xr = 1.0
    X = [S210 - S110*S300; S201 - S101*S300]
    D1 = [1-S101^2  S101*S110-S011; 
          S101*S110-S011  1-S110^2]
    
    U = sqrt(D1)  # Matrix square root
    V = U * X
    L = max(0.0, dot(V, V))
    dD1 = max(0.0, det(D1))
    R = H200 * dD1
    
    if R <= 0 || dot(X, X) < 1000*eps()
        xr = 0.0
    elseif L > R
        xr = sqrt(R / L)
    end
    
    Vr = beta * xr * V
    Xr = U \ Vr
    S210r = Xr[1] + S110*S300
    S201r = Xr[2] + S101*S300
    
    return S210r, S201r
end
