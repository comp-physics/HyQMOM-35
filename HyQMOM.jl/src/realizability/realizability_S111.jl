"""
    realizability_S111(S110, S101, S011, S210, S201, S120, S021, S102, S012, S111)

Check and correct realizability of S111 moment.

# Arguments
- Second-order standardized moments: S110, S101, S011
- Third-order standardized moments: S210, S201, S120, S021, S102, S012
- S111: Current value of S111

# Returns
- S111r: Corrected value of S111

# Algorithm
Computes bounds based on three different projections (A110, A101, A011)
and clamps S111 to the feasible range [Rmin, Rmax].
"""
@fastmath function realizability_S111(S110, S101, S011, S210, S201, S120, S021, S102, S012, S111)
    S111r = S111
    
    # Compute three alternative bounds
    A110 = ((S101-S011*S110)*S210 + (S011-S101*S110)*S120) / (1-S110^2)
    A101 = ((S110-S011*S101)*S201 + (S011-S110*S101)*S102) / (1-S101^2)
    A011 = ((S110-S101*S011)*S021 + (S101-S110*S011)*S012) / (1-S011^2)
    
    Rmin = min(A110, A101, A011)
    Rmax = max(A110, A101, A011)
    
    if S111 > Rmax
        S111r = Rmax
    elseif S111 < Rmin
        S111r = Rmin
    end
    
    return S111r
end
