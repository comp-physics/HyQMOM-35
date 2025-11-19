using LinearAlgebra: eigvals

"""
    eigenvalues6z_hyperbolic_3D(M::Vector{Float64}, flag2D::Int, Ma::Float64)

Eigenvalues of 3-D flux Jacobian in z direction.

# Arguments
- `M`: 35-element moment vector
- `flag2D`: 2D simulation flag
- `Ma`: Mach number

# Returns
- `v6min`: Minimum eigenvalue for hyperbolicity
- `v6max`: Maximum eigenvalue for hyperbolicity
- `Mr`: Corrected moment vector

# Note
M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,
     M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,
     M031,M012,M112,M013,M022]
"""
function eigenvalues6z_hyperbolic_3D(M::Vector{Float64}, flag2D::Int, Ma::Float64)
    Mr = copy(M)
    
    # CRITICAL SAFEGUARD: Check for negative density (unphysical)
    if POSITIVITY_ENABLED[] && M[1] <= 0.0
        println("  Negative/zero density detected in eigenvalues6z!")
        println("  m000 = $(M[1])")
        println("  Clamping to minimum positive value...")
        Mr[1] = max(M[1], 1e-14)
    end
    
    # Extract moments needed for WU and WV planes
    m000 = Mr[1]
    m100 = M[2]
    m200 = M[3]
    m300 = M[4]
    m400 = M[5]
    m010 = M[6]
    m020 = M[10]
    m030 = M[13]
    m040 = M[15]
    m001 = M[16]
    m101 = M[17]
    m201 = M[18]
    m301 = M[19]
    m002 = M[20]
    m102 = M[21]
    m202 = M[22]
    m003 = M[23]
    m103 = M[24]
    m004 = M[25]
    m011 = M[26]
    m021 = M[29]
    m031 = M[31]
    m012 = M[32]
    m013 = M[34]
    m022 = M[35]
    
    # SAFEGUARD: Check for negative variance before Jacobian computation
    # Compute quick estimate of variance from second moments
    if POSITIVITY_ENABLED[] && m000 > 1e-15
        var_x_est = m200/m000 - (m100/m000)^2
        var_y_est = m020/m000 - (m010/m000)^2
        var_z_est = m002/m000 - (m001/m000)^2
        
        if var_x_est < 0.0 || var_y_est < 0.0 || var_z_est < 0.0
            println("⚠️  WARNING: Negative variance detected before z-eigenvalue calculation")
            println("  var_x=$(var_x_est), var_y=$(var_y_est), var_z=$(var_z_est)")
            println("  This indicates moment realizability violation from advection")
            println("  Will attempt realizability correction...")
            
            # Apply full realizability correction
            _, _, _, Mr = Flux_closure35_and_realizable_3D(M, flag2D, Ma)
            
            # Update moments for eigenvalue calculation
            m000 = Mr[1]
            m100 = Mr[2]
            m200 = Mr[3]
            m300 = Mr[4]
            m400 = Mr[5]
            m010 = Mr[6]
            m020 = Mr[10]
            m030 = Mr[13]
            m040 = Mr[15]
            m001 = Mr[16]
            m101 = Mr[17]
            m201 = Mr[18]
            m301 = Mr[19]
            m002 = Mr[20]
            m102 = Mr[21]
            m202 = Mr[22]
            m003 = Mr[23]
            m103 = Mr[24]
            m004 = Mr[25]
            m011 = Mr[26]
            m021 = Mr[29]
            m031 = Mr[31]
            m012 = Mr[32]
            m013 = Mr[34]
            m022 = Mr[35]
        end
    end
    
    # WU moments (z-x plane)
    J6 = jacobian6(m000, m100, m200, m300, m400, m001, m101, m201, m301, m002, m102, m202, m003, m103, m004)
    lam6a = eigvals(J6)
    lam6ar = sort(real.(lam6a))
    v6min_wu = lam6ar[1]
    v6max_wu = lam6ar[6]
    v6min = v6min_wu
    v6max = v6max_wu
    
    # WV moments (z-y plane)
    # Correct mapping from MATLAB: z is first direction, y is second
    J6 = jacobian6(m000, m001, m002, m003, m004,
                   m010, m011, m012, m013,
                   m020, m021, m022,
                   m030, m031, m040)
    lam6b = eigvals(J6)
    lam6br = sort(real.(lam6b))
    v6min_wv = lam6br[1]
    v6max_wv = lam6br[6]
    v6min = min(v6min, v6min_wv)
    v6max = max(v6max, v6max_wv)
    
    # Debug: Check for pathological eigenvalues
    if abs(v6max) > 1000.0 || abs(v6min) > 1000.0
        println("WARNING in eigenvalues6z_hyperbolic_3D:")
        println("  WU (z-x): min=$(v6min_wu), max=$(v6max_wu)")
        println("  WV (z-y): min=$(v6min_wv), max=$(v6max_wv)")
        println("  Physical vel: w=$(m001/m000)")
        println("  m013 = $(m013), m031 = $(m031)")
        println("  Ratio m031/m013 = $(m031/m013)")
        println("  WV Jacobian using m031 (14th arg) = $(m031)")
        
        # Detailed diagnostics for WV plane (z-y) which seems problematic
        println("\n  Full WV moment set:")
        println("    m000=$(m000), m001=$(m001), m002=$(m002), m003=$(m003), m004=$(m004)")
        println("    m010=$(m010), m011=$(m011), m012=$(m012), m013=$(m013)")
        println("    m020=$(m020), m021=$(m021), m022=$(m022)")
        println("    m030=$(m030), m031=$(m031), m040=$(m040)")
        
        # Check for NaN/Inf in Jacobian inputs
        wv_moments = [m000, m001, m002, m003, m004, m010, m011, m012, m013, m020, m021, m022, m030, m031, m040]
        if any(!isfinite, wv_moments)
            println("  ⚠️  NON-FINITE MOMENTS IN WV PLANE!")
            for (i, val) in enumerate(wv_moments)
                if !isfinite(val)
                    println("    moment[$i] = $val")
                end
            end
        end
        
        # Check central moments
        C4, _ = M2CS4_35(M)
        C002 = C4[20]
        C020 = C4[10]
        C011 = C4[26]
        println("\n  Central moments: C002=$(C002), C020=$(C020), C011=$(C011)")
        println("  Temperature check: C002/m000=$(C002/m000), C020/m000=$(C020/m000)")
        
        println("  Will attempt correction...")
    end
    
    # Check for complex eigenvalues in z direction and correct
    if maximum(abs.(imag.(lam6a))) > 1000*eps() || maximum(abs.(imag.(lam6b))) > 1000*eps()
        # Compute mean velocities
        M000 = M[1]
        M100 = M[2]
        M010 = M[6]
        M001 = M[16]
        umean = M100/M000
        vmean = M010/M000
        wmean = M001/M000
        
        # Compute central and standardized moments
        C4, S4 = M2CS4_35(M)
        
        # Extract central moments
        C200 = C4[3]
        C300 = C4[4]
        C400 = C4[5]
        C110 = C4[7]
        C210 = C4[8]
        C310 = C4[9]
        C020 = C4[10]
        C120 = C4[11]
        C220 = C4[12]
        C030 = C4[13]
        C130 = C4[14]
        C040 = C4[15]
        C101 = C4[17]
        C201 = C4[18]
        C301 = C4[19]
        C002 = C4[20]
        C102 = C4[21]
        C202 = C4[22]
        C003 = C4[23]
        C103 = C4[24]
        C004 = C4[25]
        C011 = C4[26]
        C111 = C4[27]
        C211 = C4[28]
        C021 = C4[29]
        C121 = C4[30]
        C031 = C4[31]
        C012 = C4[32]
        C112 = C4[33]
        C013 = C4[34]
        C022 = C4[35]
        
        # Extract standardized moments
        S101 = S4[17]
        S300 = S4[4]
        S003 = S4[23]
        S400 = S4[5]
        S202 = S4[22]
        S004 = S4[25]
        S011 = S4[26]
        S030 = S4[13]
        S022 = S4[35]
        S040 = S4[15]
        
        # Force real eigenvalues for WU plane
        if maximum(abs.(imag.(lam6a))) > 1000*eps()
            S102 = S101*S003
            S201 = S101*S300
            S103 = S101*S004
            S301 = S101*S400
            s22min = (2 + 5*S003*S101*S300) / 6
            if S202 < s22min
                S202 = s22min
            end
            # Central moments from corrected standardized moments
            sC200 = sqrt(max(eps(), C200))
            sC002 = sqrt(max(eps(), C002))
            
            C201 = S201*sC200^2*sC002
            C102 = S102*sC200*sC002^2
            C301 = S301*sC200^3*sC002
            C103 = S103*sC200*sC002^3
            C202 = S202*sC200^2*sC002^2
        end
        
        # Force real eigenvalues for WV plane
        if maximum(abs.(imag.(lam6b))) > 1000*eps()
            S012 = S011*S003
            S021 = S011*S030
            S013 = S011*S004
            S031 = S011*S040
            s22min = (2 + 5*S003*S011*S030) / 6
            if S022 < s22min
                S022 = s22min
            end
            # Central moments from corrected standardized moments
            sC020 = sqrt(max(eps(), C020))
            sC002 = sqrt(max(eps(), C002))
            
            C021 = S021*sC020^2*sC002
            C012 = S012*sC020*sC002^2
            C031 = S031*sC020^3*sC002
            C013 = S013*sC020*sC002^3
            C022 = S022*sC020^2*sC002^2
        end
        
        # Integer moments from corrected central moments
        M4 = C4toM4_3D(M000, umean, vmean, wmean, C200, C110, C101, C020, C011, C002, C300,
                       C210, C201, C120, C111, C102, C030, C021, C012, C003, C400, C310, C301,
                       C220, C211, C202, C130, C121, C112, C103, C040, C031, C022, C013, C004)
        
        # Extract corrected moments
        M000 = M4[1,1,1]
        M100 = M4[2,1,1]
        M010 = M4[1,2,1]
        M001 = M4[1,1,2]
        M200 = M4[3,1,1]
        M110 = M4[2,2,1]
        M101 = M4[2,1,2]
        M020 = M4[1,3,1]
        M011 = M4[1,2,2]
        M002 = M4[1,1,3]
        M300 = M4[4,1,1]
        M210 = M4[3,2,1]
        M201 = M4[3,1,2]
        M120 = M4[2,3,1]
        M111 = M4[2,2,2]
        M102 = M4[2,1,3]
        M030 = M4[1,4,1]
        M021 = M4[1,3,2]
        M012 = M4[1,2,3]
        M003 = M4[1,1,4]
        M400 = M4[5,1,1]
        M310 = M4[4,2,1]
        M301 = M4[4,1,2]
        M220 = M4[3,3,1]
        M211 = M4[3,2,2]
        M202 = M4[3,1,3]
        M130 = M4[2,4,1]
        M121 = M4[2,3,2]
        M112 = M4[2,2,3]
        M103 = M4[2,1,4]
        M040 = M4[1,5,1]
        M031 = M4[1,4,2]
        M022 = M4[1,3,3]
        M013 = M4[1,2,4]
        M004 = M4[1,1,5]
        
        # Hyperbolic moments
        Mh = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,
              M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,
              M031,M012,M112,M013,M022]
        
        # Apply realizability constraints
        _, _, _, Mr = Flux_closure35_and_realizable_3D(Mh, flag2D, Ma)
        
        # Recompute eigenvalues with corrected moments
        m000 = Mr[1]
        m100 = Mr[2]
        m200 = Mr[3]
        m300 = Mr[4]
        m400 = Mr[5]
        m010 = Mr[6]
        m020 = Mr[10]
        m030 = Mr[13]
        m040 = Mr[15]
        m001 = Mr[16]
        m101 = Mr[17]
        m201 = Mr[18]
        m301 = Mr[19]
        m002 = Mr[20]
        m102 = Mr[21]
        m202 = Mr[22]
        m003 = Mr[23]
        m103 = Mr[24]
        m004 = Mr[25]
        m011 = Mr[26]
        m021 = Mr[29]
        m031 = Mr[31]
        m012 = Mr[32]
        m013 = Mr[34]
        m022 = Mr[35]
        
        # WU moments
        J6 = jacobian6(m000, m100, m200, m300, m400, m001, m101, m201, m301, m002, m102, m202, m003, m103, m004)
        lam6a = eigvals(J6)
        lam6ar = sort(real.(lam6a))
        v6min = lam6ar[1]
        v6max = lam6ar[6]
        
        # WV moments (z-y plane)
        # Correct mapping from MATLAB: z is first direction, y is second
        J6 = jacobian6(m000, m001, m002, m003, m004,
                       m010, m011, m012, m013,
                       m020, m021, m022,
                       m030, m031, m040)
        lam6b = eigvals(J6)
        lam6br = sort(real.(lam6b))
        v6min = min(v6min, lam6br[1])
        v6max = max(v6max, lam6br[6])
    end
    
    return v6min, v6max, Mr
end

