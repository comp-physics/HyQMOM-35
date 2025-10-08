"""
Utilities for extracting moment arrays to individual variables.

These functions replace repetitive moment unpacking code.
"""

"""
    M4_to_vars(M4)

Extract M4 array (35 moments) to individual variables.

# Returns
35 individual moment variables in canonical order.
"""
function M4_to_vars(M4::AbstractVector)
    @assert length(M4) == 35 "M4 must have 35 elements"
    
    # Order: M000, M100, M200, M300, M400, M010, M110, M210, M310, M020, M120, M220, M030, M130, M040,
    #        M001, M101, M201, M301, M002, M102, M202, M003, M103, M004, M011, M111, M211, M021, M121,
    #        M031, M012, M112, M013, M022
    return (M4[1],  M4[2],  M4[3],  M4[4],  M4[5],  M4[6],  M4[7],  M4[8],  M4[9],  M4[10],
            M4[11], M4[12], M4[13], M4[14], M4[15], M4[16], M4[17], M4[18], M4[19], M4[20],
            M4[21], M4[22], M4[23], M4[24], M4[25], M4[26], M4[27], M4[28], M4[29], M4[30],
            M4[31], M4[32], M4[33], M4[34], M4[35])
end

"""
    M4_to_vars(M4::AbstractArray{T,3})

Extract M4 3D array (5×5×5) to individual variables.

# Returns
35 individual moment variables in canonical order.
"""
function M4_to_vars(M4::AbstractArray{T,3}) where T
    @assert size(M4) == (5,5,5) "M4 must be 5×5×5"
    
    # Extract moments from 3D array using 1-based indexing
    # M[i+1, j+1, k+1] corresponds to moment M_ijk
    return (M4[1,1,1],  M4[2,1,1],  M4[3,1,1],  M4[4,1,1],  M4[5,1,1],  # M000, M100, M200, M300, M400
            M4[1,2,1],  M4[2,2,1],  M4[3,2,1],  M4[4,2,1],  # M010, M110, M210, M310
            M4[1,3,1],  M4[2,3,1],  M4[3,3,1],  # M020, M120, M220
            M4[1,4,1],  M4[2,4,1],  # M030, M130
            M4[1,5,1],  # M040
            M4[1,1,2],  M4[2,1,2],  M4[3,1,2],  M4[4,1,2],  # M001, M101, M201, M301
            M4[1,1,3],  M4[2,1,3],  M4[3,1,3],  # M002, M102, M202
            M4[1,1,4],  M4[2,1,4],  # M003, M103
            M4[1,1,5],  # M004
            M4[1,2,2],  M4[2,2,2],  M4[3,2,2],  # M011, M111, M211
            M4[1,3,2],  M4[2,3,2],  # M021, M121
            M4[1,4,2],  # M031
            M4[1,2,3],  M4[2,2,3],  # M012, M112
            M4[1,2,4],  # M013
            M4[1,3,3])  # M022
end

"""
    M5_to_vars(M5)

Extract M5 array (56 moments) to individual variables.

# Returns
56 individual moment variables in canonical order (up to 5th order).
"""
function M5_to_vars(M5::AbstractVector)
    @assert length(M5) == 56 "M5 must have 56 elements"
    
    # 35 moments from M4 + 21 fifth-order moments
    return (M5[1],  M5[2],  M5[3],  M5[4],  M5[5],  M5[6],  M5[7],  M5[8],  M5[9],  M5[10],
            M5[11], M5[12], M5[13], M5[14], M5[15], M5[16], M5[17], M5[18], M5[19], M5[20],
            M5[21], M5[22], M5[23], M5[24], M5[25], M5[26], M5[27], M5[28], M5[29], M5[30],
            M5[31], M5[32], M5[33], M5[34], M5[35], M5[36], M5[37], M5[38], M5[39], M5[40],
            M5[41], M5[42], M5[43], M5[44], M5[45], M5[46], M5[47], M5[48], M5[49], M5[50],
            M5[51], M5[52], M5[53], M5[54], M5[55], M5[56])
end

"""
    M5_to_vars(M5::AbstractArray{T,3})

Extract M5 3D array (6×6×6) to individual variables.

# Returns
56 individual moment variables in canonical order (up to 5th order).
"""
function M5_to_vars(M5::AbstractArray{T,3}) where T
    @assert size(M5) == (6,6,6) "M5 must be 6×6×6"
    
    # Extract moments from 3D array using 1-based indexing
    # M[i+1, j+1, k+1] corresponds to moment M_ijk
    # Order matches MATLAB extract_M5_array exactly
    
    # Order 0-1
    M000 = M5[1,1,1]
    M100 = M5[2,1,1]; M010 = M5[1,2,1]; M001 = M5[1,1,2]
    
    # Order 2
    M200 = M5[3,1,1]; M110 = M5[2,2,1]; M101 = M5[2,1,2]
    M020 = M5[1,3,1]; M011 = M5[1,2,2]; M002 = M5[1,1,3]
    
    # Order 3
    M300 = M5[4,1,1]; M210 = M5[3,2,1]; M201 = M5[3,1,2]
    M120 = M5[2,3,1]; M111 = M5[2,2,2]; M102 = M5[2,1,3]
    M030 = M5[1,4,1]; M021 = M5[1,3,2]; M012 = M5[1,2,3]; M003 = M5[1,1,4]
    
    # Order 4
    M400 = M5[5,1,1]; M310 = M5[4,2,1]; M301 = M5[4,1,2]
    M220 = M5[3,3,1]; M211 = M5[3,2,2]; M202 = M5[3,1,3]
    M130 = M5[2,4,1]; M121 = M5[2,3,2]; M112 = M5[2,2,3]; M103 = M5[2,1,4]
    M040 = M5[1,5,1]; M031 = M5[1,4,2]; M022 = M5[1,3,3]; M013 = M5[1,2,4]; M004 = M5[1,1,5]
    
    # Order 5 (closure)
    M500 = M5[6,1,1]
    M410 = M5[5,2,1]; M320 = M5[4,3,1]; M230 = M5[3,4,1]; M140 = M5[2,5,1]
    M401 = M5[5,1,2]; M302 = M5[4,1,3]; M203 = M5[3,1,4]; M104 = M5[2,1,5]
    M311 = M5[4,2,2]; M221 = M5[3,3,2]; M131 = M5[2,4,2]
    M212 = M5[3,2,3]; M113 = M5[2,2,4]; M122 = M5[2,3,3]
    M050 = M5[1,6,1]; M041 = M5[1,5,2]; M032 = M5[1,4,3]; M023 = M5[1,3,4]; M014 = M5[1,2,5]; M005 = M5[1,1,6]
    
    return (M000, M100, M010, M001, M200, M110, M101, M020, M011, M002,
            M300, M210, M201, M120, M111, M102, M030, M021, M012, M003,
            M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004,
            M500, M410, M320, M230, M140, M401, M302, M203, M104, M311, M221, M131, M212, M113, M122,
            M050, M041, M032, M023, M014, M005)
end

"""
    S_to_C_batch(S110, S101, S011, S300, ..., sC200, sC020, sC002)

Batch convert standardized moments to central moments.

# Algorithm
For each moment: `C_ijk = S_ijk * sC200^i * sC020^j * sC002^k`

# Returns
All central moments in the same order as inputs.
"""
function S_to_C_batch(S110, S101, S011, S300, S210, S201, S120, S111, S102, S030, S021, S012, S003,
                      S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004,
                      S500, S410, S401, S320, S311, S302, S230, S221, S212, S203, S140, S131, S122, S113, S104,
                      S050, S041, S032, S023, S014, S005,
                      sC200, sC020, sC002)
    
    # 2nd order
    C110 = S110 * sC200 * sC020
    C101 = S101 * sC200 * sC002
    C011 = S011 * sC020 * sC002
    
    # 3rd order
    C300 = S300 * sC200^3
    C210 = S210 * sC200^2 * sC020
    C201 = S201 * sC200^2 * sC002
    C120 = S120 * sC200 * sC020^2
    C111 = S111 * sC200 * sC020 * sC002
    C102 = S102 * sC200 * sC002^2
    C030 = S030 * sC020^3
    C021 = S021 * sC020^2 * sC002
    C012 = S012 * sC020 * sC002^2
    C003 = S003 * sC002^3
    
    # 4th order
    C400 = S400 * sC200^4
    C310 = S310 * sC200^3 * sC020
    C301 = S301 * sC200^3 * sC002
    C220 = S220 * sC200^2 * sC020^2
    C211 = S211 * sC200^2 * sC020 * sC002
    C202 = S202 * sC200^2 * sC002^2
    C130 = S130 * sC200 * sC020^3
    C121 = S121 * sC200 * sC020^2 * sC002
    C112 = S112 * sC200 * sC020 * sC002^2
    C103 = S103 * sC200 * sC002^3
    C040 = S040 * sC020^4
    C031 = S031 * sC020^3 * sC002
    C022 = S022 * sC020^2 * sC002^2
    C013 = S013 * sC020 * sC002^3
    C004 = S004 * sC002^4
    
    # 5th order
    C500 = S500 * sC200^5
    C410 = S410 * sC200^4 * sC020
    C401 = S401 * sC200^4 * sC002
    C320 = S320 * sC200^3 * sC020^2
    C311 = S311 * sC200^3 * sC020 * sC002
    C302 = S302 * sC200^3 * sC002^2
    C230 = S230 * sC200^2 * sC020^3
    C221 = S221 * sC200^2 * sC020^2 * sC002
    C212 = S212 * sC200^2 * sC020 * sC002^2
    C203 = S203 * sC200^2 * sC002^3
    C140 = S140 * sC200 * sC020^4
    C131 = S131 * sC200 * sC020^3 * sC002
    C122 = S122 * sC200 * sC020^2 * sC002^2
    C113 = S113 * sC200 * sC020 * sC002^3
    C104 = S104 * sC200 * sC002^4
    C050 = S050 * sC020^5
    C041 = S041 * sC020^4 * sC002
    C032 = S032 * sC020^3 * sC002^2
    C023 = S023 * sC020^2 * sC002^3
    C014 = S014 * sC020 * sC002^4
    C005 = S005 * sC002^5
    
    return (C110, C101, C011, C300, C210, C201, C120, C111, C102, C030, C021, C012, C003,
            C400, C310, C301, C220, C211, C202, C130, C121, C112, C103, C040, C031, C022, C013, C004,
            C500, C410, C401, C320, C311, C302, C230, C221, C212, C203, C140, C131, C122, C113, C104,
            C050, C041, C032, C023, C014, C005)
end
