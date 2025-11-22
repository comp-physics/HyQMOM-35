"""
Diagnostic utilities for moment checking and validation.
"""

"""
    check2D(S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)

Check and correct 2D moments in 3D code.

# Algorithm
If `|S11| >= 1`, collapses to boundary case.
Otherwise, applies `realizability(:2D, ...)` to correct cross moments.

# Returns
10 corrected moments
"""
function check2D(S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)
    h2min = 1000 * eps()
    
    if abs(S11) >= 1
        # Collapse both S11 >= 1 and S11 <= -1 branches
        S11 = sign(S11)
        s3m = sqrt(abs(S30 * S03))
        s4m = sqrt(S40 * S04)
        H2m = max(h2min, s4m - s3m^2 - 1)
        s4m = H2m + s3m^2 + 1
        S30 = sign(S30) * s3m
        S03 = S11 * S30
        S40 = s4m
        S04 = s4m
        S12 = S11 * S03
        S21 = S11 * S30
        S13 = S11 * S04
        S31 = S11 * S40
        S22 = s4m
    else
        # Check and correct realizability of cross moments
        S21, S12, S31, S22, S13 = realizability(Symbol("2D"), S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)
    end
    
    return (S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)
end

"""
    check2D_all_planes(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                       S101, S201, S301, S102, S202, S003, S103, S004,
                       S011, S021, S031, S012, S022, S013)

Apply check2D realizability to all three coordinate planes.

# Returns
30 corrected moments (10 per plane: XY, XZ, YZ)

# Algorithm
Applies `check2D` to each of the three coordinate planes independently.
"""
function check2D_all_planes(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                            S101, S201, S301, S102, S202, S003, S103, S004,
                            S011, S021, S031, S012, S022, S013)
    
    # Apply check2D to XY plane
    S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1 =
        check2D(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040)
    
    # Apply check2D to XZ plane
    S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2 =
        check2D(S300, S400, S101, S201, S301, S102, S202, S003, S103, S004)
    
    # Apply check2D to YZ plane
    S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3 =
        check2D(S030, S040, S011, S021, S031, S012, S022, S003, S013, S004)
    
    return (S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1,
            S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2,
            S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3)
end