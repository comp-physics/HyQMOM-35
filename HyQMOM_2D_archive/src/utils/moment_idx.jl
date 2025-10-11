"""
    moment_idx(order)

Return linear indices for moment arrays.

Helper function for accessing moments in 3D arrays.
"""
function moment_idx(order::Int)
    if order == 4
        # 35 moments up to 4th order
        # Manually computed sub2ind([5 5 5], i, j, k) = i + 5*(j-1) + 25*(k-1)
        return [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 16, 17, 21,
                26, 27, 28, 29, 51, 52, 53, 76, 77, 101, 31, 32, 33, 36, 37,
                41, 56, 57, 81, 61]
    elseif order == 5
        # 56 moments up to 5th order
        error("5th order indices not yet implemented")
    else
        error("Order must be 4 or 5")
    end
end

"""
    moment_idx(name)

Return the index of a moment by name (e.g., "M110" -> 7).

# Arguments
- `name`: Moment name as a string (e.g., "M110", "M200")

# Returns
- Integer index in the 35-element moment vector
"""
function moment_idx(name::String)
    # Parse moment name (e.g., "M110" -> i=1, j=1, k=0)
    if length(name) != 4 || name[1] != 'M'
        error("Invalid moment name: $name. Expected format: Mijk")
    end
    
    i = parse(Int, name[2:2]) + 1  # Convert to 1-based
    j = parse(Int, name[3:3]) + 1
    k = parse(Int, name[4:4]) + 1
    
    # Compute linear index in 5x5x5 array
    linear_idx = i + 5*(j-1) + 25*(k-1)
    
    # Map to position in 35-element vector
    idx_c = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 16, 17, 21,
             26, 27, 28, 29, 51, 52, 53, 76, 77, 101, 31, 32, 33, 36, 37,
             41, 56, 57, 81, 61]
    
    pos = findfirst(==(linear_idx), idx_c)
    if pos === nothing
        error("Moment $name is not in the 35-element set")
    end
    
    return pos
end

"""
    moment_names(order)

Return canonical moment names for given order.

# Arguments
- `order`: Moment order (4 or 5)

# Returns
- Vector of moment name strings
"""
function moment_names(order::Int)
    if order == 4
        return ["M000", "M100", "M200", "M300", "M400",
                "M010", "M110", "M210", "M310", "M020", "M120", "M220", "M030", "M130", "M040",
                "M001", "M101", "M201", "M301", "M002", "M102", "M202", "M003", "M103", "M004",
                "M011", "M111", "M211", "M021", "M121", "M031", "M012", "M112", "M013", "M022"]
    elseif order == 5
        return ["M000", "M100", "M200", "M300", "M400", "M500",
                "M010", "M110", "M210", "M310", "M410",
                "M020", "M120", "M220", "M320",
                "M030", "M130", "M230",
                "M040", "M140",
                "M050",
                "M001", "M101", "M201", "M301", "M401",
                "M011", "M111", "M211", "M311",
                "M021", "M121", "M221",
                "M031", "M131",
                "M041",
                "M002", "M102", "M202", "M302",
                "M012", "M112", "M212",
                "M022", "M122",
                "M032",
                "M003", "M103", "M203",
                "M013", "M113",
                "M023",
                "M004", "M104",
                "M014",
                "M005"]
    else
        error("Order must be 4 or 5")
    end
end
