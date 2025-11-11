"""
Initial Condition Utilities for HyQMOM

Provides flexible setup for initial moment distributions including:
- Uniform background states
- Cubic/box regions with specified velocities and densities
- Coordinate-based placement (physical or index-based)
"""

"""
    CubicRegion

Defines a cubic/box region with specified moment properties.

# Fields
- `center::NTuple{3,Float64}`: Center position (x, y, z) in physical coordinates
- `width::NTuple{3,Float64}`: Width in each direction (dx, dy, dz)
- `density::Float64`: Density (rho) in the region
- `velocity::NTuple{3,Float64}`: Mean velocity (U, V, W)
- `temperature::Float64`: Temperature (for covariance matrix)

# Example
```julia
# Cube centered at origin, 0.2 units wide, moving in +x direction
cube = CubicRegion(
    center = (0.0, 0.0, 0.0),
    width = (0.2, 0.2, 0.2),
    density = 1.0,
    velocity = (0.5, 0.0, 0.0),
    temperature = 1.0
)
```
"""
struct CubicRegion
    center::NTuple{3,Float64}
    width::NTuple{3,Float64}
    density::Float64
    velocity::NTuple{3,Float64}
    temperature::Float64
end

"""
    CubicRegion(; center, width, density, velocity, temperature=1.0)

Keyword constructor for CubicRegion.
"""
function CubicRegion(; center, width, density, velocity, temperature=1.0)
    return CubicRegion(
        Tuple(Float64.(center)),
        Tuple(Float64.(width)),
        Float64(density),
        Tuple(Float64.(velocity)),
        Float64(temperature)
    )
end

"""
    initialize_moment_field(grid_params, background, regions; 
                           r110=0.0, r101=0.0, r011=0.0)

Initialize a 3D moment field with background state and cubic regions.

# Arguments
- `grid_params`: NamedTuple with (Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax, xm, ym, zm)
- `background`: CubicRegion defining the background state (uniform)
- `regions`: Vector of CubicRegion objects to place in the domain
- `r110`, `r101`, `r011`: Correlation coefficients for covariance matrix

# Returns
- `M`: 4D array (nx, ny, nz, Nmom) of initialized moments

# Example
```julia
# Background at rest with low density
background = CubicRegion(
    center = (0.0, 0.0, 0.0),
    width = (Inf, Inf, Inf),  # Fills entire domain
    density = 0.01,
    velocity = (0.0, 0.0, 0.0),
    temperature = 1.0
)

# Two crossing jets
jet1 = CubicRegion(
    center = (-0.25, -0.25, 0.0),
    width = (0.1, 0.1, 0.1),
    density = 1.0,
    velocity = (0.707, 0.707, 0.0),  # Moving northeast
    temperature = 1.0
)

jet2 = CubicRegion(
    center = (0.25, 0.25, 0.0),
    width = (0.1, 0.1, 0.1),
    density = 1.0,
    velocity = (-0.707, -0.707, 0.0),  # Moving southwest
    temperature = 1.0
)

M = initialize_moment_field(grid_params, background, [jet1, jet2])
```
"""
function initialize_moment_field(grid_params, background::CubicRegion, 
                                regions::Vector{CubicRegion};
                                r110=0.0, r101=0.0, r011=0.0)
    
    Nx = grid_params.Nx
    Ny = grid_params.Ny
    Nz = grid_params.Nz
    xm = grid_params.xm
    ym = grid_params.ym
    zm = grid_params.zm
    Nmom = 35
    
    M = zeros(Float64, Nx, Ny, Nz, Nmom)
    
    # Build covariance matrix (same for all regions, using temperature)
    # Each region can have different temperature
    
    # Initialize with background
    for k in 1:Nz
        for j in 1:Ny
            for i in 1:Nx
                M[i, j, k, :] = region_to_moments(background, r110, r101, r011)
            end
        end
    end
    
    # Place each cubic region
    for region in regions
        place_cubic_region!(M, region, xm, ym, zm, r110, r101, r011)
    end
    
    return M
end

"""
    initialize_moment_field_mpi(decomp, grid_params, background, regions; kwargs...)

MPI-aware version that only initializes local subdomain.

# Arguments
- `decomp`: MPI decomposition structure
- `grid_params`: Grid parameters
- `background`: Background CubicRegion
- `regions`: Vector of CubicRegion objects
- `halo`: Halo size (default: 2)

# Returns
- `M`: 4D array with halos (nx+2*halo, ny+2*halo, nz, Nmom)
"""
function initialize_moment_field_mpi(decomp, grid_params, background::CubicRegion,
                                    regions::Vector{CubicRegion};
                                    r110=0.0, r101=0.0, r011=0.0, halo=2)
    
    Nx = grid_params.Nx
    Ny = grid_params.Ny
    Nz = grid_params.Nz
    xm = grid_params.xm
    ym = grid_params.ym
    zm = grid_params.zm
    Nmom = 35
    
    # Local grid size
    nx = decomp.local_size[1]
    ny = decomp.local_size[2]
    nz = decomp.local_size[3]
    
    # Global index ranges
    i0i1 = decomp.istart_iend
    j0j1 = decomp.jstart_jend
    k0k1 = decomp.kstart_kend
    
    # Allocate with halos
    M = zeros(Float64, nx + 2*halo, ny + 2*halo, nz, Nmom)
    
    # Background moments
    Mr_bg = region_to_moments(background, r110, r101, r011)
    
    # Calculate grid spacing for overlap detection
    dx = Nx > 1 ? abs(xm[2] - xm[1]) : 0.0
    dy = Ny > 1 ? abs(ym[2] - ym[1]) : 0.0
    dz = Nz > 1 ? abs(zm[2] - zm[1]) : 0.0
    cell_size = (dx, dy, dz)
    
    # Fill local subdomain
    for kk in 1:nz
        gk = k0k1[1] + kk - 1  # global k index
        z_coord = zm[gk]
        
        for ii in 1:nx
            gi = i0i1[1] + ii - 1  # global i index
            x_coord = xm[gi]
            
            for jj in 1:ny
                gj = j0j1[1] + jj - 1  # global j index
                y_coord = ym[gj]
                
                # Default: background
                Mr = Mr_bg
                
                # Check each region (later regions overwrite earlier ones)
                # Use cell overlap detection for robustness with narrow jets
                cell_center = (x_coord, y_coord, z_coord)
                for region in regions
                    if cell_overlaps_cube(cell_center, cell_size, region)
                        Mr = region_to_moments(region, r110, r101, r011)
                    end
                end
                
                M[ii + halo, jj + halo, kk, :] = Mr
            end
        end
    end
    
    return M
end

"""
    region_to_moments(region::CubicRegion, r110, r101, r011)

Convert a CubicRegion specification to a moment vector.
"""
function region_to_moments(region::CubicRegion, r110, r101, r011)
    U, V, W = region.velocity
    rho = region.density
    T = region.temperature
    
    # Covariance matrix components
    C200 = T
    C020 = T
    C002 = T
    C110 = r110 * sqrt(C200 * C020)
    C101 = r101 * sqrt(C200 * C002)
    C011 = r011 * sqrt(C020 * C002)
    
    return InitializeM4_35(rho, U, V, W, C200, C110, C101, C020, C011, C002)
end

"""
    point_in_cube(point::NTuple{3,Float64}, cube::CubicRegion)

Check if a point (x, y, z) is inside a cubic region.

Note: This checks if the **point** (typically a grid cell center) is inside
the cube. For narrow jets on coarse grids, jets may be missed if their width
is smaller than the grid spacing. To avoid this, ensure jet width >= grid spacing
or use cell_overlaps_cube() for more robust overlap detection.
"""
function point_in_cube(point::NTuple{3,Float64}, cube::CubicRegion)
    x, y, z = point
    cx, cy, cz = cube.center
    wx, wy, wz = cube.width
    
    # Handle infinite width (background)
    wx = isinf(wx) ? Inf : wx / 2.0
    wy = isinf(wy) ? Inf : wy / 2.0
    wz = isinf(wz) ? Inf : wz / 2.0
    
    return (abs(x - cx) <= wx &&
            abs(y - cy) <= wy &&
            abs(z - cz) <= wz)
end

"""
    cell_overlaps_cube(cell_center::NTuple{3,Float64}, cell_size::NTuple{3,Float64}, cube::CubicRegion)

Check if a grid cell (defined by center and size) overlaps with a cubic region.

More robust than point_in_cube() for detecting when narrow jets intersect coarse grid cells.
A cell and cube overlap if they intersect in all three dimensions.

# Arguments
- `cell_center`: (x, y, z) coordinates of cell center
- `cell_size`: (dx, dy, dz) cell dimensions
- `cube`: CubicRegion to check overlap with

# Returns
- `true` if cell and cube overlap in all three dimensions
"""
function cell_overlaps_cube(cell_center::NTuple{3,Float64}, cell_size::NTuple{3,Float64}, cube::CubicRegion)
    x, y, z = cell_center
    dx, dy, dz = cell_size
    cx, cy, cz = cube.center
    wx, wy, wz = cube.width
    
    # Handle infinite width (background always overlaps)
    if isinf(wx) || isinf(wy) || isinf(wz)
        return true
    end
    
    # Check overlap in each dimension independently
    # Two 1D intervals [a, b] and [c, d] overlap if: max(a,c) < min(b,d)
    # Cell interval: [center - size/2, center + size/2]
    # Cube interval: [center - width/2, center + width/2]
    
    x_overlap = (max(x - dx/2, cx - wx/2) < min(x + dx/2, cx + wx/2))
    y_overlap = (max(y - dy/2, cy - wy/2) < min(y + dy/2, cy + wy/2))
    z_overlap = (max(z - dz/2, cz - wz/2) < min(z + dz/2, cz + wz/2))
    
    return x_overlap && y_overlap && z_overlap
end

"""
    place_cubic_region!(M, region::CubicRegion, xm, ym, zm, r110, r101, r011; use_overlap=true)

Place a cubic region into an existing moment field (modifies M in-place).

# Arguments
- `M`: Moment field array (Nx, Ny, Nz, Nmom)
- `region`: CubicRegion to place
- `xm`, `ym`, `zm`: Grid cell centers
- `r110`, `r101`, `r011`: Correlation coefficients
- `use_overlap`: If true, use cell overlap detection (more robust). If false, use point-in-cube (legacy)
"""
function place_cubic_region!(M, region::CubicRegion, xm, ym, zm, r110, r101, r011; use_overlap=true)
    Nx = length(xm)
    Ny = length(ym)
    Nz = length(zm)
    
    Mr = region_to_moments(region, r110, r101, r011)
    
    if use_overlap
        # Calculate grid spacing (assume uniform spacing)
        dx = Nx > 1 ? abs(xm[2] - xm[1]) : 0.0
        dy = Ny > 1 ? abs(ym[2] - ym[1]) : 0.0
        dz = Nz > 1 ? abs(zm[2] - zm[1]) : 0.0
        
        for k in 1:Nz
            for j in 1:Ny
                for i in 1:Nx
                    cell_center = (xm[i], ym[j], zm[k])
                    cell_size = (dx, dy, dz)
                    if cell_overlaps_cube(cell_center, cell_size, region)
                        M[i, j, k, :] = Mr
                    end
                end
            end
        end
    else
        # Legacy: point-in-cube check (may miss narrow jets on coarse grids)
        for k in 1:Nz
            for j in 1:Ny
                for i in 1:Nx
                    point = (xm[i], ym[j], zm[k])
                    if point_in_cube(point, region)
                        M[i, j, k, :] = Mr
                    end
                end
            end
        end
    end
end

"""
    crossing_jets_ic(Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax;
                    Ma=0.0, rhol=1.0, rhor=0.01, T=1.0, 
                    jet_size=0.1, jet_offset=0.0)

Create standard crossing jets initial condition using the flexible system.

This is a convenience function that creates the classic two-jet configuration.

# Arguments
- Grid parameters: `Nx`, `Ny`, `Nz`, domain bounds
- `Ma`: Mach number (determines jet velocity magnitude)
- `rhol`: Jet density
- `rhor`: Background density
- `T`: Temperature
- `jet_size`: Fraction of domain for jet width (default: 0.1 = 10%)
- `jet_offset`: Offset from center in units of jet_size (default: 0)

# Returns
- `background::CubicRegion`
- `regions::Vector{CubicRegion}` with two jets
"""
function crossing_jets_ic(Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax;
                         Ma=0.0, rhol=1.0, rhor=0.01, T=1.0, 
                         jet_size=0.1, jet_offset=0.0)
    
    # Domain center and size
    x_center = (xmin + xmax) / 2.0
    y_center = (ymin + ymax) / 2.0
    z_center = (zmin + zmax) / 2.0
    
    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin
    
    # Jet dimensions
    jet_width_x = jet_size * Lx
    jet_width_y = jet_size * Ly
    jet_width_z = jet_size * Lz
    
    # Jet positions (offset from center)
    offset_dist = jet_offset * jet_width_x
    
    # Crossing jets velocity
    Uc = Ma / sqrt(2.0)
    
    # Background at rest
    background = CubicRegion(
        center = (x_center, y_center, z_center),
        width = (Inf, Inf, Inf),
        density = rhor,
        velocity = (0.0, 0.0, 0.0),
        temperature = T
    )
    
    # Bottom-left jet (moving northeast)
    jet1 = CubicRegion(
        center = (x_center - offset_dist - jet_width_x/2, 
                  y_center - offset_dist - jet_width_y/2, 
                  z_center),
        width = (jet_width_x, jet_width_y, jet_width_z),
        density = rhol,
        velocity = (Uc, Uc, 0.0),
        temperature = T
    )
    
    # Top-right jet (moving southwest)
    jet2 = CubicRegion(
        center = (x_center + offset_dist + jet_width_x/2,
                  y_center + offset_dist + jet_width_y/2,
                  z_center),
        width = (jet_width_x, jet_width_y, jet_width_z),
        density = rhol,
        velocity = (-Uc, -Uc, 0.0),
        temperature = T
    )
    
    return background, [jet1, jet2]
end

