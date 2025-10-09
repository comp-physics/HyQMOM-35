"""
    realizability_S310(S110, S101, S011, S300, S210, S201, S120, S111, S310, S220, H200, beta)

Check and correct realizability of S310 and S220 moments.

# Arguments
- S110, S101, S011: Second-order standardized moments
- S300, S210, S201, S120, S111: Third-order moments
- S310, S220: Fourth-order moments to correct
- H200: Variance-related quantity
- beta: Scaling factor (typically 1.0)

# Returns
- S310r, S220r: Corrected moments

# Algorithm
Computes bounds for S310 based on complex quadratic forms involving
both 2nd and 3rd order moments. The bound depends on S220 and H200.
"""
@fastmath function realizability_S310(S110, S101, S011, S300, S210, S201, S120, S111, S310, S220, H200, beta)
    S310r = S310
    S220r = S220
    
    # Compute b310 coefficient
    b310 = S111*((1 - S110^2)*S201 + (S101*S110 - S011)*S210 + (S011*S110 - S101)*S300) +
           S120*((1 - S101^2)*S210 + (S011*S101 - S110)*S300 + (S101*S110 - S011)*S201) +
           S210*((1 - S011^2)*S300 + (S011*S101 - S110)*S210 + (S011*S110 - S101)*S201)
    
    # Build D1 matrix and compute L1
    D1 = [1-S101^2  S101*S110-S011; 
          S101*S110-S011  1-S110^2]
    dD1 = det(D1)
    V1 = [S210 - S110*S300; S201 - S101*S300]
    L1 = dot(V1, D1 * V1) / dD1
    
    # Compute G310b
    G310b = H200 - L1
    if G310b < 0
        G310b = 0.0
    end
    
    # Build D2 matrix and compute L2
    D2 = [1-S011^2      S011*S101-S110  S011*S110-S101;
          S011*S101-S110  1-S101^2      S101*S110-S011;
          S011*S110-S101  S101*S110-S011  1-S110^2]
    V2 = [S210; S120; S111]
    L2 = dot(V2, D2 * V2) / dD1
    
    # Compute G310a
    G310a = S220 - S110^2 - L2 + 1000*eps()
    if G310a < 0 || L2 < 0
        G310a = 0.0
    end
    
    # Compute bounds
    G310 = G310a * G310b
    sG310 = beta * sqrt(G310)
    
    s310min = S110 + b310/dD1 - sG310
    s310max = S110 + b310/dD1 + sG310
    
    if S310 <= s310min
        S310r = s310min
    elseif S310 >= s310max
        S310r = s310max
    end
    
    return S310r, S220r
end
