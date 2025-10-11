"""
    realizable_2D(S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)

Find realizable moment set in 2D.

# Arguments
- S30, S40: Third and fourth-order moments along first axis
- S11: Second-order cross-moment
- S21, S31: Third-order cross-moments
- S12, S22: Cross-moments
- S03, S13, S04: Moments along second axis

# Returns
- S21, S12, S31, S22, S13: Corrected realizable moments

# Algorithm
Sequentially enforces realizability constraints:
1. Check realizability of S11 (first minor)
2. Check realizability of S12 and S21
3. Check realizability of S22 (second minor)
4. Check realizability of S13 and S31 given S22
5. Check third minor using Cholesky determinant and polynomial roots
"""
function realizable_2D(S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)
    # Check realizability of S11
    Del1 = max(0.0, 1 - S11^2)
    
    H20 = max(eps(), S40 - S30^2 - 1)
    H02 = max(eps(), S04 - S03^2 - 1)
    
    # Check realizability of S12 and S21
    G1 = sqrt(Del1 * H02)
    s12min = S11*S03 - G1
    s12max = S11*S03 + G1
    if S12 <= s12min
        S12 = s12min
    elseif S12 >= s12max
        S12 = s12max
    end
    
    G1 = sqrt(Del1 * H20)
    s21min = S11*S30 - G1
    s21max = S11*S30 + G1
    if S21 <= s21min
        S21 = s21min
    elseif S21 >= s21max
        S21 = s21max
    end
    
    # At this point first minor is nonnegative
    
    # Check realizability of S22
    G22 = sqrt((H20 + S30^2) * (H02 + S03^2))
    s22min = max(S11^2, 1 - G22)
    s22max = 1 + G22
    if S22 < s22min
        S22 = s22min
    elseif S22 > s22max
        S22 = s22max
    end
    
    # Given S22, check realizability of S13 and S31
    G31 = (Del1*S22 - Del1*S11^2 - S12^2 + 2*S11*S12*S21 - S21^2) * 
          (Del1*H20 - (S21 - S11*S30)^2)
    G13 = (Del1*S22 - Del1*S11^2 - S12^2 + 2*S11*S12*S21 - S21^2) * 
          (Del1*H02 - (S12 - S11*S03)^2)
    
    if G31 < 0 || G13 < 0
        G31 = 0.0
        G13 = 0.0
    else
        G13 = sqrt(G13)
        G31 = sqrt(G31)
    end
    
    s31min = S11 + (S12*S21 + S21*S30 - S11*S21^2 - S11*S12*S30 - G31) / (Del1 + eps())
    s31max = S11 + (S12*S21 + S21*S30 - S11*S21^2 - S11*S12*S30 + G31) / (Del1 + eps())
    if S31 <= s31min
        S31 = s31min
    elseif S31 >= s31max
        S31 = s31max
    end
    
    s13min = S11 + (S12*S21 + S12*S03 - S11*S12^2 - S11*S21*S03 - G13) / (Del1 + eps())
    s13max = S11 + (S12*S21 + S12*S03 - S11*S12^2 - S11*S21*S03 + G13) / (Del1 + eps())
    if S13 <= s13min
        S13 = s13min
    elseif S13 >= s13max
        S13 = s13max
    end
    
    # At this point second minor is nonnegative, check third for S22
    L3 = delta2starchol_L3(S03, S04, S11, S12, S13, S21, S22, S30, S31, S40)
    if L3 < 0
        # Check realizability of S22 using roots of degree 3 polynomial
        R = rootsR(Del1, H02, H20, S03, S11, S12, S13, S21, S30, S31)
        Rr = sort(real(R))
        s22maxr = min(s22max, Rr[3])
        s22minr = max(s22min, Rr[2])
        if s22maxr < s22minr
            s22minr = s22maxr
        end
        S22a = S22
        if S22 > s22maxr
            S22a = s22maxr
        elseif S22 < s22minr
            S22a = s22minr
        end
        # Check for complex roots, in which case L3 remains < 0
        if maximum(abs.(imag(R))) / maximum(abs.(R)) > 1000*eps()
            # Complex roots: do not change S22
        else
            S22 = S22a
        end
    end
    
    return S21, S12, S31, S22, S13
end
