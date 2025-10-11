"""
    Flux_closure35_and_realizable_3D(M4, flag2D, Ma)

Compute 3D fluxes for all moments and correct unrealizable moments.

This is the main closure function that orchestrates the entire moment pipeline:
1. Convert raw moments to central and standardized moments
2. Enforce univariate realizability
3. Check and correct 2D realizability in all planes
4. Handle edge/corner cases
5. Apply 3D realizability
6. Compute 5th-order closure via HyQMOM
7. Convert back to raw moments
8. Assemble flux vectors

# Arguments
- `M4`: 35-element moment vector up to 4th order
  `[M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,
    M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,
    M031,M012,M112,M013,M022]`
- `flag2D`: 2D simulation flag (1 for 2D, 0 for 3D)
- `Ma`: Mach number

# Returns
- `Fx`: X-direction flux moments (35 elements, 5th order)
- `Fy`: Y-direction flux moments (35 elements, 5th order)
- `Fz`: Z-direction flux moments (35 elements, 5th order)
- `M4r`: Realizable 4th-order moments (35 elements)

# Algorithm
See Flux_closure35_and_realizable_3D.m for detailed algorithm description.
"""
function Flux_closure35_and_realizable_3D(M4::AbstractVector, flag2D::Int, Ma::Real; debug_label="", debug_output=false)
    # Constants
    s3max = 4.0 + abs(Ma) / 2.0
    h2min = 1.0e-8
    itrealmax = 6
    
    # Extract basic quantities
    M000 = M4[1]
    umean = M4[2] / M000   # M100/M000
    vmean = M4[6] / M000   # M010/M000
    wmean = M4[16] / M000  # M001/M000
    
    # DEBUG: Check if M[3] is around 0.04 (the problematic input)
    debug_this = debug_output && abs(M4[3] - 0.04) < 0.01 && M4[1] > 0.03 && M4[1] < 0.04
    if debug_this
        println("\n[DEBUG] Flux_closure35_and_realizable_3D $(debug_label):")
        @printf("  Input M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", 
                M4[1], M4[2], M4[3], M4[4], M4[5])
    end
    
    # Compute central and standardized moments
    C4, S4 = M2CS4_35(M4)
    
    # Extract key central moments (variances)
    C200 = max(eps(), C4[3])
    C020 = max(eps(), C4[10])
    C002 = max(eps(), C4[20])
    
    # Extract standardized moments
    S300=S4[4];  S400=S4[5];  S110=S4[7];  S210=S4[8];  S310=S4[9]
    S120=S4[11]; S220=S4[12]; S030=S4[13]; S130=S4[14]; S040=S4[15]
    S101=S4[17]; S201=S4[18]; S301=S4[19]; S102=S4[21]; S202=S4[22]
    S003=S4[23]; S103=S4[24]; S004=S4[25]; S011=S4[26]; S111=S4[27]
    S211=S4[28]; S021=S4[29]; S121=S4[30]; S031=S4[31]; S012=S4[32]
    S112=S4[33]; S013=S4[34]; S022=S4[35]
    
    ## Check univariate moments
    S300, S400, H200 = enforce_univariate(S300, S400, h2min, s3max)
    S030, S040, H020 = enforce_univariate(S030, S040, h2min, s3max)
    S003, S004, H002 = enforce_univariate(S003, S004, h2min, s3max)
    
    ## 4th-order moments: check maximum bounds on S220, S202, S022
    A220 = sqrt((H200 + S300^2) * (H020 + S030^2))
    S220 = realizability(:S220, S110, S220, A220)
    A202 = sqrt((H200 + S300^2) * (H002 + S003^2))
    S202 = realizability(:S220, S101, S202, A202)
    A022 = sqrt((H020 + S030^2) * (H002 + S003^2))
    S022 = realizability(:S220, S011, S022, A022)
    
    ## Check and correct realizability of 2D moments
    (S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1,
     S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2,
     S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3) =
        check2D_all_planes(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                          S101, S201, S301, S102, S202, S003, S103, S004,
                          S011, S021, S031, S012, S022, S013)
    
    # Store original values
    S111r = S111
    S211r = S211
    S121r = S121
    S112r = S112
    S300r = S300
    S030r = S030
    S003r = S003
    S400r = S400
    S040r = S040
    S004r = S004
    
    ## Check for non-realizable 2D moments outside the box
    R110 = 1 - S110^2
    R101 = 1 - S101^2
    R011 = 1 - S011^2
    
    if R110 <= 0 || R101 <= 0 || R011 <= 0
        # Treat cases where one or more 2D 2nd-order moments is non-realizable (corners and edges)
        (S110, S101, S011, S300, S030, S003, S400, S040, S004,
         S210, S201, S120, S111, S102, S021, S012,
         S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013) =
            edge_corner_correction(R110, R101, R011,
                                  S110r, S101r, S011r,
                                  S300r1, S030r1, S400r1, S040r1,
                                  S300r2, S003r2, S400r2, S004r2,
                                  S030r3, S003r3, S040r3, S004r3,
                                  S300r, S030r, S003r, S400r, S040r, S004r,
                                  S210r, S201r, S120r, S111r, S102r, S021r, S012r,
                                  S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r,
                                  S031r, S022r, S013r)
        # Recheck 2D cross moments
        S210, S120, S310, S220, S130 = realizability(Symbol("2D"), S300, S400, S110, S210, S310, S120, S220, S030, S130, S040)
        S201, S301, S102, S202, S103 = realizability(Symbol("2D"), S300, S400, S101, S201, S301, S102, S202, S003, S103, S004)
        S021, S031, S012, S022, S013 = realizability(Symbol("2D"), S030, S040, S011, S021, S031, S012, S022, S003, S013, S004)
    else
        # Treat cases with faces or interior of 2nd-order moment space
        (S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
         S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
         S211, S021, S121, S031, S012, S112, S013, S022, flag220) =
            realizability(Symbol("3D"), S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                         S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
                         S211, S021, S121, S031, S012, S112, S013, S022)
        itreal = 0
        while flag220 == 1 && itreal < itrealmax
            itreal += 1
            # Repeat if S220 has changed
            (S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
             S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
             S211, S021, S121, S031, S012, S112, S013, S022, _) =
                realizability(Symbol("3D"), S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                             S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
                             S211, S021, S121, S031, S012, S112, S013, S022)
        end
    end
    
    ## Force moments for pure 2D case with S011=S101=0
    if flag2D == 1
        if abs(S011) > h2min || abs(S101) > h2min
            @warn "flag2D==1 but S011 or S101 is nonzero!"
        end
        S011 = 0.0
        S101 = 0.0
        S201 = 0.0
        S021 = 0.0
        S012 = 0.0
        S102 = 0.0
        S111 = 0.0
        S202 = 1.0
        S022 = 1.0
        S112 = S110
        S121 = 0.0
        S211 = 0.0
        S103 = 0.0
        S013 = 0.0
        S301 = 0.0
        S031 = 0.0
    end
    
    ## 3D HyQMOM closures for 5th-order standardized moments
    (S500, S410, S320, S230, S140, S401, S302, S203, S104, S311,
     S221, S131, S212, S113, S122, S050, S041, S032, S023, S014, S005) =
        hyqmom_3D(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                 S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
                 S211, S021, S121, S031, S012, S112, S013, S022)
    
    ## 5th-order central moments from corrected standardized moments
    sC200 = sqrt(max(eps(), C200))
    sC020 = sqrt(max(eps(), C020))
    sC002 = sqrt(max(eps(), C002))
    
    # Batch convert S->C
    (C110, C101, C011, C300, C210, C201, C120, C111, C102, C030, C021, C012, C003,
     C400, C310, C301, C220, C211, C202, C130, C121, C112, C103, C040, C031, C022, C013, C004,
     C500, C410, C401, C320, C311, C302, C230, C221, C212, C203, C140, C131, C122, C113, C104,
     C050, C041, C032, C023, C014, C005) =
        S_to_C_batch(S110, S101, S011, S300, S210, S201, S120, S111, S102, S030, S021, S012, S003,
                    S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004,
                    S500, S410, S401, S320, S311, S302, S230, S221, S212, S203, S140, S131, S122, S113, S104,
                    S050, S041, S032, S023, S014, S005,
                    sC200, sC020, sC002)
    
    ## 5th-order moments from central moments
    M5 = C5toM5_3D(M000, umean, vmean, wmean, C200, C110, C101, C020, C011, C002,
                   C300, C210, C201, C120, C111, C102, C030, C021, C012, C003,
                   C400, C310, C301, C220, C211, C202, C130, C121, C112, C103, C040, C031, C022, C013, C004,
                   C500, C410, C320, C230, C140, C401, C302, C203, C104, C311, C221, C131, C212, C113, C122,
                   C050, C041, C032, C023, C014, C005)
    
    # Extract M5 array to individual variables
    (M000, M100, M010, M001, M200, M110, M101, M020, M011, M002,
     M300, M210, M201, M120, M111, M102, M030, M021, M012, M003,
     M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004,
     M500, M410, M320, M230, M140, M401, M302, M203, M104, M311, M221, M131, M212, M113, M122,
     M050, M041, M032, M023, M014, M005) = M5_to_vars(M5)
    
    ## Flux closures
    Fx = [M100,M200,M300,M400,M500,M110,M210,M310,M410,M120,M220,M320,M130,M230,M140,
          M101,M201,M301,M401,M102,M202,M302,M103,M203,M104,M111,M211,M311,M121,M221,
          M131,M112,M212,M113,M122]
    Fy = [M010,M110,M210,M310,M410,M020,M120,M220,M320,M030,M130,M230,M040,M140,M050,
          M011,M111,M211,M311,M012,M112,M212,M013,M113,M014,M021,M121,M221,M031,M131,
          M041,M022,M122,M023,M032]
    Fz = [M001,M101,M201,M301,M401,M011,M111,M211,M311,M021,M121,M221,M031,M131,M041,
          M002,M102,M202,M302,M003,M103,M203,M004,M104,M005,M012,M112,M212,M022,M122,
          M032,M013,M113,M014,M023]
    
    # Realizable moments
    M4r = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,
           M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,
           M031,M012,M112,M013,M022]
    
    # DEBUG: Check output
    if debug_this
        @printf("  Output M4r[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", 
                M4r[1], M4r[2], M4r[3], M4r[4], M4r[5])
        @printf("  M[3] change: %.6e -> %.6e (%.2fx)\n", M4[3], M4r[3], M4r[3] / M4[3])
        @printf("  M4r element 3 is M200 = %.6e\n", M200)
    end
    
    # Make a defensive copy to avoid aliasing issues
    M4r_copy = copy(M4r)
    
    return Fx, Fy, Fz, M4r_copy
end
