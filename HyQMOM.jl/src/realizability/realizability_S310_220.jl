"""
    realizability_S310_220(S110, S101, S011, S210, S120, S111, S220)

Check and correct realizability of S220 for S310 constraint.

# Arguments
- S110, S101, S011: Second-order standardized moments
- S210, S120, S111: Third-order moments
- S220: Fourth-order moment to correct

# Returns
- S220r: Corrected S220 (lower bound)

# Algorithm
Computes a lower bound for S220 based on the constraint that
certain quadratic forms must be non-negative.
"""
@fastmath function realizability_S310_220(S110, S101, S011, S210, S120, S111, S220)
    S220r = S220
    
    D1 = [1-S101^2  S101*S110-S011; 
          S101*S110-S011  1-S110^2]
    dD1 = det(D1)
    
    D2 = [1-S011^2      S011*S101-S110  S011*S110-S101;
          S011*S101-S110  1-S101^2      S101*S110-S011;
          S011*S110-S101  S101*S110-S011  1-S110^2]
    V2 = [S210; S120; S111]
    L2 = dot(V2, D2 * V2) / dD1
    S220_min = S110^2 + L2 + 1000*eps()
    
    if S220_min > S220r
        S220r = S220_min
    end
    
    return S220r
end
