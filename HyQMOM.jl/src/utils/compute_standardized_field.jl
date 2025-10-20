"""
    compute_standardized_field(M)

Compute standardized moments for entire field.

# Arguments
- `M`: 4D array (Nx, Ny, Nz, Nmom) of raw moments

# Returns
- `S`: 4D array (Nx, Ny, Nz, Nmom) of standardized moments

# Notes
Applies M2CS4_35 to each grid point to extract standardized moments.
This is useful for visualization and analysis of snapshot data.
"""
function compute_standardized_field(M::Array{Float64,4})
    Nx, Ny, Nz, Nmom = size(M)
    @assert Nmom == 35 "Expected 35 moments, got $Nmom"
    
    S = zeros(Float64, Nx, Ny, Nz, Nmom)
    
    for k in 1:Nz
        for j in 1:Ny
            for i in 1:Nx
                mom = M[i, j, k, :]
                _, S4 = M2CS4_35(mom)
                S[i, j, k, :] = S4
            end
        end
    end
    
    return S
end

"""
    compute_central_field(M)

Compute central moments for entire field.

# Arguments
- `M`: 4D array (Nx, Ny, Nz, Nmom) of raw moments

# Returns
- `C`: 4D array (Nx, Ny, Nz, Nmom) of central moments

# Notes
Applies M2CS4_35 to each grid point to extract central moments.
"""
function compute_central_field(M::Array{Float64,4})
    Nx, Ny, Nz, Nmom = size(M)
    @assert Nmom == 35 "Expected 35 moments, got $Nmom"
    
    C = zeros(Float64, Nx, Ny, Nz, Nmom)
    
    for k in 1:Nz
        for j in 1:Ny
            for i in 1:Nx
                mom = M[i, j, k, :]
                C4, _ = M2CS4_35(mom)
                C[i, j, k, :] = C4
            end
        end
    end
    
    return C
end

"""
    get_standardized_moment(S, moment_name::String)

Extract a specific standardized moment from the field.

# Arguments
- `S`: 4D array (Nx, Ny, Nz, Nmom) of standardized moments
- `moment_name`: Name of moment (e.g., "S110", "S101", "S022")

# Returns
- 3D array (Nx, Ny, Nz) of the requested moment

# Supported moments
Second order: S110, S101, S011
Third order: S300, S210, S201, S120, S111, S102, S030, S021, S012, S003
Fourth order: S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004
"""
function get_standardized_moment(S::Array{Float64,4}, moment_name::String)
    # Moment name to index mapping (matches the 35-element vector ordering)
    moment_map = Dict(
        "S000" => 1,
        "S100" => 2, "S200" => 3, "S300" => 4, "S400" => 5,
        "S010" => 6, "S110" => 7, "S210" => 8, "S310" => 9,
        "S020" => 10, "S120" => 11, "S220" => 12,
        "S030" => 13, "S130" => 14,
        "S040" => 15,
        "S001" => 16, "S101" => 17, "S201" => 18, "S301" => 19,
        "S002" => 20, "S102" => 21, "S202" => 22,
        "S003" => 23, "S103" => 24,
        "S004" => 25,
        "S011" => 26, "S111" => 27, "S211" => 28,
        "S021" => 29, "S121" => 30,
        "S031" => 31,
        "S012" => 32, "S112" => 33,
        "S013" => 34,
        "S022" => 35
    )
    
    if !haskey(moment_map, moment_name)
        error("Unknown moment name: $moment_name. Use format like 'S110', 'S101', etc.")
    end
    
    idx = moment_map[moment_name]
    return S[:, :, :, idx]
end

"""
    get_central_moment(C, moment_name::String)

Extract a specific central moment from the field.

# Arguments
- `C`: 4D array (Nx, Ny, Nz, Nmom) of central moments
- `moment_name`: Name of moment (e.g., "C110", "C101", "C022")

# Returns
- 3D array (Nx, Ny, Nz) of the requested moment

# Supported moments
Same indexing as standardized moments but with 'C' prefix.
"""
function get_central_moment(C::Array{Float64,4}, moment_name::String)
    # Convert C to S in the name and use the same mapping
    S_name = replace(moment_name, "C" => "S")
    return get_standardized_moment(C, S_name)  # Uses same indexing
end

