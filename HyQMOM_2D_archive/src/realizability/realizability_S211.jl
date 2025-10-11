"""
    realizability_S211(e11, e22, e33, e12, e13, d23, S211, beta)

Check and correct realizability of S211 moment.

# Arguments
- e11, e22, e33, e12, e13: Elements from realizability matrix
- d23: Offset term
- S211: Current value of S211
- beta: Scaling factor (typically 1.0)

# Returns
- S211r: Corrected value of S211

# Algorithm
Computes bounds based on quadratic form involving e-matrix elements.
Clamps S211 to feasible range [s211min, s211max].
"""
@fastmath function realizability_S211(e11, e22, e33, e12, e13, d23, S211, beta)
    S211r = S211
    b211 = e12 * e13
    G211 = max(0.0, (e11*e22 - e13^2) * (e11*e33 - e12^2))
    sG211 = beta * sqrt(G211)
    s211min = d23 + (b211 - sG211) / e11
    s211max = d23 + (b211 + sG211) / e11
    
    if S211 <= s211min
        S211r = s211min
    elseif S211 >= s211max
        S211r = s211max
    end
    
    return S211r
end
