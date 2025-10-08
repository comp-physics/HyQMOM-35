"""
    eigenvalues6_hyperbolic_3D(M, axis, flag2D, Ma)

Compute unified eigenvalues of 3D flux Jacobian.

# Arguments
- `M`: 35-element moment vector
- `axis`: Direction (1=X direction, 2=Y direction)
- `flag2D`: 2D simulation flag
- `Ma`: Mach number

# Returns
- `v6min`, `v6max`: Min/max eigenvalues for hyperbolicity
- `Mr`: Corrected moment vector

# Algorithm
1. Compute eigenvalues from current moments
2. Check for complex eigenvalues
3. If complex, correct moments to ensure real eigenvalues
4. Recompute eigenvalues with corrected moments
"""
function eigenvalues6_hyperbolic_3D(M::AbstractVector, axis::Int, flag2D::Int, Ma::Real)
    Mr = copy(M)
    
    # DEBUG
    debug_this = abs(M[3] - 0.04) < 0.01 && M[1] > 0.03 && M[1] < 0.04
    if debug_this
        @printf("    eigenvalues6_hyperbolic_3D axis=%d: Input M[3]=%.6e\n", axis, M[3])
    end
    
    # Compute eigenvalues from current moments
    v6min, v6max, lam6a, lam6b = compute_eigenvalues_for_axis(M, axis)
    
    # Check for complex eigenvalues and correct if needed
    if maximum(abs.(imag(lam6a))) > 1000*eps() || maximum(abs.(imag(lam6b))) > 1000*eps()
        if debug_this
            println("    → Complex eigenvalues detected, correcting...")
        end
        # Correct moments to ensure real eigenvalues
        M_corrected = correct_moments_for_real_eigenvalues(M, axis, lam6a, lam6b, flag2D, Ma)
        
        if debug_this
            @printf("    → After correction: M_corrected[3]=%.6e\n", M_corrected[3])
        end
        
        # Recompute eigenvalues with corrected moments
        v6min, v6max, _, _ = compute_eigenvalues_for_axis(M_corrected, axis)
        Mr = M_corrected
    end
    
    if debug_this
        @printf("    eigenvalues6_hyperbolic_3D axis=%d: Output Mr[3]=%.6e\n", axis, Mr[3])
    end
    
    return v6min, v6max, Mr
end

"""
    compute_eigenvalues_for_axis(M, axis)

Helper: Compute eigenvalues for given axis.
"""
function compute_eigenvalues_for_axis(M::AbstractVector, axis::Int)
    # Use axis_moment_slice to extract the right moment vectors
    if axis == 1  # X direction: UV and UW planes
        moments_uv = axis_moment_slice(M, 1)  # UV plane
        moments_uw = axis_moment_slice(M, 3)  # UW plane
        return compute_jacobian_eigenvalues(moments_uv, moments_uw)
    else  # Y direction: VU and VW planes
        moments_vu = axis_moment_slice(M, 2)  # VU plane
        moments_vw = axis_moment_slice(M, 4)  # VW plane
        return compute_jacobian_eigenvalues(moments_vu, moments_vw)
    end
end

"""
    correct_moments_for_real_eigenvalues(M, axis, lam6a, lam6b, flag2D, Ma)

Helper: Correct moments to ensure real eigenvalues.
"""
function correct_moments_for_real_eigenvalues(M::AbstractVector, axis::Int, lam6a, lam6b, flag2D::Int, Ma::Real)
    # Compute mean velocities
    M000 = M[1];  M100 = M[2];  M010 = M[6];  M001 = M[16]
    umean = M100 / M000
    vmean = M010 / M000
    wmean = M001 / M000
    
    # Compute central and standardized moments
    C4, S4 = M2CS4_35(M)
    
    # Extract necessary moments from C4 and S4 vectors
    C200 = C4[3];   C020 = C4[10];  C002 = C4[20]
    C110 = C4[7];   C101 = C4[17];  C011 = C4[26]
    
    # Extract ALL C4 moments (we'll update some and need to pass all to C4toM4_3D)
    C_all = extract_all_C4_moments(C4)
    
    # Extract key standardized moments
    S110 = S4[7];   S300 = S4[4];   S030 = S4[13]
    S400 = S4[5];   S220 = S4[12];  S040 = S4[15]
    S101 = S4[17];  S202 = S4[22];  S003 = S4[23]
    S004 = S4[25];  S011 = S4[26];  S022 = S4[35]
    
    # Correct for first eigenvalue pair (lam6a) - always UV plane
    if maximum(abs.(imag(lam6a))) > 1000*eps()
        C_all, S220 = correct_uv_plane(C_all, S110, S300, S030, S400, S040, S220, C200, C020)
    end
    
    # Correct for second eigenvalue pair (lam6b) - depends on axis
    if maximum(abs.(imag(lam6b))) > 1000*eps()
        if axis == 1  # X: UW moments
            C_all = correct_uw_plane(C_all, S101, S300, S003, S400, S004, S202, C200, C002)
        else  # Y: VW moments
            C_all = correct_vw_plane(C_all, S011, S030, S003, S040, S004, S022, C020, C002)
        end
    end
    
    # Reconstruct M4 from corrected central moments
    M4 = C4toM4_3D(M000, umean, vmean, wmean,
                   C_all.C200, C_all.C110, C_all.C101, C_all.C020, C_all.C011, C_all.C002,
                   C_all.C300, C_all.C210, C_all.C201, C_all.C120, C_all.C111, C_all.C102,
                   C_all.C030, C_all.C021, C_all.C012, C_all.C003,
                   C_all.C400, C_all.C310, C_all.C301, C_all.C220, C_all.C211, C_all.C202,
                   C_all.C130, C_all.C121, C_all.C112, C_all.C103,
                   C_all.C040, C_all.C031, C_all.C022, C_all.C013, C_all.C004)
    
    # Extract to vector
    M000, M100, M010, M001, M200, M110, M101, M020, M011, M002,
    M300, M210, M201, M120, M111, M102, M030, M021, M012, M003,
    M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004 =
        M4_to_vars(M4)
    
    Mh = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,
          M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,
          M031,M012,M112,M013,M022]
    
    # Apply realizability constraints
    _, _, _, M_corrected = Flux_closure35_and_realizable_3D(Mh, flag2D, Ma, debug_label="[eigenvalues6_nested]")
    
    return M_corrected
end

"""
    extract_all_C4_moments(C4)

Helper: Extract all C4 moments into a named tuple.
"""
function extract_all_C4_moments(C4::AbstractVector)
    return (C200=C4[3],   C300=C4[4],   C400=C4[5],
            C110=C4[7],   C210=C4[8],   C310=C4[9],
            C020=C4[10],  C120=C4[11],  C220=C4[12],
            C030=C4[13],  C130=C4[14],  C040=C4[15],
            C101=C4[17],  C201=C4[18],  C301=C4[19],
            C002=C4[20],  C102=C4[21],  C202=C4[22],
            C003=C4[23],  C103=C4[24],  C004=C4[25],
            C011=C4[26],  C111=C4[27],  C211=C4[28],
            C021=C4[29],  C121=C4[30],  C031=C4[31],
            C012=C4[32],  C112=C4[33],  C013=C4[34],
            C022=C4[35])
end

"""
    correct_uv_plane(C, S110, S300, S030, S400, S040, S220, C200, C020)

Helper: Correct UV plane moments.
"""
function correct_uv_plane(C, S110, S300, S030, S400, S040, S220, C200, C020)
    S120 = S110 * S030
    S210 = S110 * S300
    S130 = S110 * S040
    S310 = S110 * S400
    s22min = (2 + 5*S030*S110*S300) / 6
    if S220 < s22min
        S220 = s22min
    end
    
    sC200 = sqrt(max(eps(), C200))
    sC020 = sqrt(max(eps(), C020))
    
    # Return new named tuple with updated values
    C_new = merge(C, (C210 = S210 * sC200^2 * sC020,
                      C120 = S120 * sC200 * sC020^2,
                      C310 = S310 * sC200^3 * sC020,
                      C130 = S130 * sC200 * sC020^3,
                      C220 = S220 * sC200^2 * sC020^2))
    
    return C_new, S220
end

"""
    correct_uw_plane(C, S101, S300, S003, S400, S004, S202, C200, C002)

Helper: Correct UW plane moments.
"""
function correct_uw_plane(C, S101, S300, S003, S400, S004, S202, C200, C002)
    S102 = S101 * S003
    S201 = S101 * S300
    S103 = S101 * S004
    S301 = S101 * S400
    s22min = (2 + 5*S003*S101*S300) / 6
    if S202 < s22min
        S202 = s22min
    end
    
    sC200 = sqrt(max(eps(), C200))
    sC002 = sqrt(max(eps(), C002))
    
    return merge(C, (C201 = S201 * sC200^2 * sC002,
                     C102 = S102 * sC200 * sC002^2,
                     C301 = S301 * sC200^3 * sC002,
                     C103 = S103 * sC200 * sC002^3,
                     C202 = S202 * sC200^2 * sC002^2))
end

"""
    correct_vw_plane(C, S011, S030, S003, S040, S004, S022, C020, C002)

Helper: Correct VW plane moments.
"""
function correct_vw_plane(C, S011, S030, S003, S040, S004, S022, C020, C002)
    S012 = S011 * S003
    S021 = S011 * S030
    S013 = S011 * S004
    S031 = S011 * S040
    s22min = (2 + 5*S003*S011*S030) / 6
    if S022 < s22min
        S022 = s22min
    end
    
    sC020 = sqrt(max(eps(), C020))
    sC002 = sqrt(max(eps(), C002))
    
    return merge(C, (C021 = S021 * sC020^2 * sC002,
                     C012 = S012 * sC020 * sC002^2,
                     C031 = S031 * sC020^3 * sC002,
                     C013 = S013 * sC020 * sC002^3,
                     C022 = S022 * sC020^2 * sC002^2))
end
