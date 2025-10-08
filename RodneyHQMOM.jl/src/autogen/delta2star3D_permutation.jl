"""
    delta2star3D_permutation(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                             S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
                             S211, S021, S121, S031, S012, S112, S013, S022)

Find 6 permutations of Δ₂* (all must be positive for realizability).

This function computes the Δ₂* matrix for all 6 permutations of the coordinate axes:
- E1: Δ₂*_ijk (original)
- E2: Δ₂*_ikj (swap j and k)
- E3: Δ₂*_jik (swap i and j)
- E4: Δ₂*_jki (cyclic permutation)
- E5: Δ₂*_kij (cyclic permutation)
- E6: Δ₂*_kji (swap i and k)

# Arguments
All 28 standardized moments up to 4th order (excluding those determined by symmetry).

# Returns
- `E1, E2, E3, E4, E5, E6`: Six 6×6 matrices, one for each permutation

# Realizability
For moments to be realizable, all six matrices must have non-negative determinants
for all principal minors.
"""
function delta2star3D_permutation(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                                   S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
                                   S211, S021, S121, S031, S012, S112, S013, S022)
    
    # varDelta2star_ijk (original ordering)
    E1 = delta2star3D(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                      S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
                      S211, S021, S121, S031, S012, S112, S013, S022)
    
    # varDelta2star_ikj (swap j and k)
    E2 = delta2star3D(S300, S400, S101, S201, S301, S102, S202, S003, S103, S004,
                      S110, S210, S310, S120, S220, S030, S130, S040, S011, S111,
                      S211, S012, S112, S013, S021, S121, S031, S022)
    
    # varDelta2star_jik (swap i and j)
    E3 = delta2star3D(S030, S040, S110, S120, S130, S210, S220, S300, S310, S400,
                      S011, S021, S031, S012, S022, S003, S013, S004, S101, S111,
                      S121, S201, S211, S301, S102, S112, S103, S202)
    
    # varDelta2star_jki (cyclic permutation: i→j, j→k, k→i)
    E4 = delta2star3D(S003, S004, S101, S102, S103, S201, S202, S300, S301, S400,
                      S011, S012, S013, S021, S022, S030, S031, S040, S110, S111,
                      S112, S210, S211, S310, S120, S121, S130, S220)
    
    # varDelta2star_kij (cyclic permutation: i→k, k→j, j→i)
    E5 = delta2star3D(S030, S040, S011, S021, S031, S012, S022, S003, S013, S004,
                      S110, S120, S130, S210, S220, S300, S310, S400, S101, S111,
                      S121, S102, S112, S103, S201, S211, S301, S202)
    
    # varDelta2star_kji (swap i and k)
    E6 = delta2star3D(S003, S004, S011, S012, S013, S021, S022, S030, S031, S040,
                      S101, S102, S103, S201, S202, S300, S301, S400, S110, S111,
                      S112, S120, S121, S130, S210, S211, S310, S220)
    
    return E1, E2, E3, E4, E5, E6
end
