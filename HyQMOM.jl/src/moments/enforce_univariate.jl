"""
    enforce_univariate(S3, S4, h2min, s3max)

Enforce realizability bounds on univariate moments.

Ensures H = S4 - S3² - 1 ≥ h2min and |S3| ≤ s3max.

# Arguments
- `S3`: Third standardized moment
- `S4`: Fourth standardized moment
- `h2min`: Minimum allowed H value
- `s3max`: Maximum allowed |S3| value

# Returns
- `S3`: Corrected third moment
- `S4`: Corrected fourth moment
- `H`: Corrected H value
"""
function enforce_univariate(S3, S4, h2min, s3max)
    # Compute H
    H = S4 - S3^2 - 1
    
    # Enforce minimum H
    if H <= h2min
        H = h2min
        S4 = H + S3^2 + 1
    end
    
    # Enforce S3 bounds
    if S3 < -s3max
        S3 = -s3max
        S4 = H + S3^2 + 1
    elseif S3 > s3max
        S3 = s3max
        S4 = H + S3^2 + 1
    end
    
    return S3, S4, H
end
