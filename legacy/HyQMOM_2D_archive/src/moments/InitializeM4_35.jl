"""
    InitializeM4_35(M000, umean, vmean, wmean, C200, C110, C101, C020, C011, C002)

Compute 3D fourth-order joint Gaussian moments from physical parameters.

# Arguments
- `M000`: Number density
- `umean, vmean, wmean`: Mean velocities (M100/M000, M010/M000, M001/M000)
- `C200, C110, C101, C020, C011, C002`: Covariance matrix elements

# Returns
- `M`: 35-element vector of raw moments up to 4th order
"""
function InitializeM4_35(M000, umean, vmean, wmean, C200, C110, C101, C020, C011, C002)
    # Standardized moments for Maxwellian (Gaussian)
    # 3rd order: all zero (Gaussian is symmetric)
    S300=0.0; S210=0.0; S201=0.0; S120=0.0; S111=0.0; S102=0.0
    S030=0.0; S021=0.0; S012=0.0; S003=0.0
    
    # 4th order: diagonal = 3 (Gaussian kurtosis), cross = 1 (independent variables)
    S400=3.0; S310=0.0; S301=0.0; S220=1.0; S211=0.0; S202=1.0
    S130=0.0; S121=0.0; S112=0.0; S103=0.0; S040=3.0; S031=0.0
    S022=1.0; S013=0.0; S004=3.0
    
    # Compute central moments from standardized moments
    C4 = S4toC4_3D_r(C200, C110, C101, C020, C011, C002,
                     S300, S210, S201, S120, S111, S102, S030, S021, S012, S003,
                     S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004)
    
    # Extract central moments from 3D array (order matches M4_to_vars for 3D arrays)
    C000, C100, C200, C300, C400,
    C010, C110, C210, C310,
    C020, C120, C220,
    C030, C130,
    C040,
    C001, C101, C201, C301,
    C002, C102, C202,
    C003, C103,
    C004,
    C011, C111, C211,
    C021, C121,
    C031,
    C012, C112,
    C013,
    C022 = M4_to_vars(C4)
    
    # Compute raw moments from central moments
    M4 = C4toM4_3D(M000, umean, vmean, wmean,
                   C200, C110, C101, C020, C011, C002,
                   C300, C210, C201, C120, C111, C102, C030, C021, C012, C003,
                   C400, C310, C301, C220, C211, C202, C130, C121, C112, C103, C040, C031, C022, C013, C004)
    
    # Extract raw moments from 3D array (order matches M4_to_vars for 3D arrays)
    M000, M100, M200, M300, M400,
    M010, M110, M210, M310,
    M020, M120, M220,
    M030, M130,
    M040,
    M001, M101, M201, M301,
    M002, M102, M202,
    M003, M103,
    M004,
    M011, M111, M211,
    M021, M121,
    M031,
    M012, M112,
    M013,
    M022 = M4_to_vars(M4)
    
    # Pack into 35-element vector (standard ordering)
    M = [M000, M100, M200, M300, M400, M010, M110, M210, M310, M020, M120, M220, M030, M130, M040,
         M001, M101, M201, M301, M002, M102, M202, M003, M103, M004, M011, M111, M211, M021, M121,
         M031, M012, M112, M013, M022]
    
    return M
end
