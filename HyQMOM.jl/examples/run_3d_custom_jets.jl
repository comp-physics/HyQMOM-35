"""
3D Custom Jet Configuration with Flexible Initial Conditions

Demonstrates how to create custom initial conditions with:
- Arbitrary number of cubic regions
- Custom positions (centroids)
- Custom sizes (widths)
- Custom velocities in each cube
- Fully controllable via command-line or code

Usage:
  # Use predefined configurations
  julia --project=. examples/run_3d_custom_jets.jl --config triple-jet
  julia --project=. examples/run_3d_custom_jets.jl --config quad-jet
  
  # Standard crossing jets (same as run_3d_jets_timeseries.jl)
  julia --project=. examples/run_3d_custom_jets.jl --config crossing
  
  # With MPI
  mpiexec -n 4 julia --project=. examples/run_3d_custom_jets.jl --config triple-jet --Np 100
  
  # Override physics parameters
  julia --project=. examples/run_3d_custom_jets.jl --config quad-jet --Ma 1.5 --tmax 0.3

To create your own configuration, edit the `get_jet_configuration` function below.
"""

# Load GLMakie first to enable visualization support in HyQMOM
try
    import GLMakie
catch e
    @warn """
    GLMakie is not installed. Visualization will not be available.
    To install: julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'
    """
end

using HyQMOM
using MPI
using Printf

# Load parameter parsing utilities
include("parse_params.jl")

"""
    get_jet_configuration(config_name, params)

Define different jet configurations.

Returns (background::CubicRegion, jets::Vector{CubicRegion})
"""
function get_jet_configuration(config_name::String, params)
    # Extract domain parameters
    xmin, xmax = params.xmin, params.xmax
    ymin, ymax = params.ymin, params.ymax
    zmin, zmax = params.zmin, params.zmax
    
    x_center = (xmin + xmax) / 2
    y_center = (ymin + ymax) / 2
    z_center = (zmin + zmax) / 2
    
    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin
    
    # Common parameters
    Ma = params.Ma
    rhol = 1.0  # Jet density
    rhor = 0.01  # Background density
    T = params.T
    
    # Jet size (10% of domain)
    jet_width = 0.1 * Lx
    
    # Background (uniform, low density, at rest)
    background = CubicRegion(
        center = (x_center, y_center, z_center),
        width = (Inf, Inf, Inf),
        density = rhor,
        velocity = (0.0, 0.0, 0.0),
        temperature = T
    )
    
    if config_name == "crossing"
        # Classic crossing jets (moving diagonally in x-y plane)
        Uc = Ma / sqrt(2.0)
        offset = jet_width * 0.6
        
        jets = [
            CubicRegion(
                center = (x_center - offset, y_center - offset, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (Uc, Uc, 0.0),  # Moving ↗
                temperature = T
            ),
            CubicRegion(
                center = (x_center + offset, y_center + offset, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (-Uc, -Uc, 0.0),  # Moving ↙
                temperature = T
            )
        ]

    elseif config_name == "crossing3D"
        # Classic crossing jets (moving diagonally in x-y plane)
        Uc = Ma / sqrt(2.0)
        offset = jet_width * 0.6
        
        jets = [
            CubicRegion(
                center = (x_center - offset, y_center - offset, z_center - offset),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (Uc, Uc, Uc),  # Moving ↗
                temperature = T
            ),
            CubicRegion(
                center = (x_center + offset, y_center + offset, z_center+offset),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (-Uc, -Uc, -Uc),  # Moving ↙
                temperature = T
            )
        ]
        
    elseif config_name == "triple-jet"
        # Three jets: one central, two opposing
        Uc = Ma
        offset = jet_width * 1.2
        
        jets = [
            # Central jet (moving in +x)
            CubicRegion(
                center = (x_center, y_center, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (Uc, 0.0, 0.0),
                temperature = T
            ),
            # Left jet (moving in +y)
            CubicRegion(
                center = (x_center - offset, y_center, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (0.0, Uc, 0.0),
                temperature = T
            ),
            # Right jet (moving in -y)
            CubicRegion(
                center = (x_center + offset, y_center, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (0.0, -Uc, 0.0),
                temperature = T
            )
        ]
        
    elseif config_name == "quad-jet"
        # Four jets in a square, all moving toward center
        Uc = Ma / sqrt(2.0)
        offset = jet_width * 1.2
        
        jets = [
            # Bottom-left (moving ↗)
            CubicRegion(
                center = (x_center - offset, y_center - offset, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (Uc, Uc, 0.0),
                temperature = T
            ),
            # Bottom-right (moving ↖)
            CubicRegion(
                center = (x_center + offset, y_center - offset, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (-Uc, Uc, 0.0),
                temperature = T
            ),
            # Top-left (moving ↘)
            CubicRegion(
                center = (x_center - offset, y_center + offset, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (Uc, -Uc, 0.0),
                temperature = T
            ),
            # Top-right (moving ↙)
            CubicRegion(
                center = (x_center + offset, y_center + offset, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (-Uc, -Uc, 0.0),
                temperature = T
            )
        ]
        
    elseif config_name == "vertical-jet"
        # Single jet moving in +z direction
        jets = [
            CubicRegion(
                center = (x_center, y_center, z_center - 0.2 * Lz),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (0.0, 0.0, Ma),
                temperature = T
            )
        ]
        
    elseif config_name == "spiral"
        # Three jets arranged for spiral motion
        Uc = Ma / sqrt(2.0)
        offset = jet_width * 1.0
        
        jets = [
            # Tangential velocities for rotation
            CubicRegion(
                center = (x_center + offset, y_center, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (0.0, Uc, Uc),  # +y, +z
                temperature = T
            ),
            CubicRegion(
                center = (x_center, y_center + offset, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (-Uc, 0.0, Uc),  # -x, +z
                temperature = T
            ),
            CubicRegion(
                center = (x_center, y_center, z_center),
                width = (jet_width, jet_width, jet_width),
                density = rhol,
                velocity = (Uc, Uc, 0.0),  # +x, +y
                temperature = T
            )
        ]
        
    else
        error("Unknown configuration: $config_name. Options: crossing, triple-jet, quad-jet, vertical-jet, spiral")
    end
    
    return background, jets
end

# Main execution
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Parse parameters with custom config option
params = parse_simulation_params(
    # Code defaults (can be overridden with command-line args)
    Np = 20,
    Nz = 20,
    tmax = 0.1,
    Ma = 0.0,
    Kn = 1.0,
    CFL = 0.3,
    xmin = -0.5,
    xmax = 0.5,
    ymin = -0.5,
    ymax = 0.5,
    zmin = -0.5,
    zmax = 0.5,
    snapshot_interval = 1,
    homogeneous_z = false,
    config = "crossing"  # Default configuration
)

# Print parameter summary
if rank == 0
    println("="^70)
    println("3D CUSTOM JET CONFIGURATION")
    println("="^70)
end

print_params_summary(params, rank=rank, comm=comm)

# Get jet configuration
if rank == 0
    println("\nJet Configuration: $(params.config)")
    println("="^70)
end

# Get jet configuration
background, jets = get_jet_configuration(params.config, params)

# Show configuration details
if rank == 0
    println("\nBackground:")
    println("  Density: $(background.density)")
    println("  Velocity: $(background.velocity)")
    
    println("\nJets: $(length(jets))")
    for (i, jet) in enumerate(jets)
        println("  Jet $i:")
        println("    Center: $(jet.center)")
        println("    Width: $(jet.width)")
        println("    Density: $(jet.density)")
        println("    Velocity: $(jet.velocity)")
    end
    println("="^70)
end

# Store IC configuration in params for simulation_runner to use
params_with_ic = merge(params, (
    ic_background = background,
    ic_jets = jets,
    use_custom_ic = true
))

# Run simulation with custom initial conditions
if params.snapshot_interval > 0
    # With snapshots
    if rank == 0
        println("\nRunning with snapshot collection...")
    end
    
    snapshots, grid = run_simulation_with_snapshots(params_with_ic; 
                                                     snapshot_interval=params.snapshot_interval)
    
    if rank == 0 && snapshots !== nothing
        println("\n" * "="^70)
        println("SIMULATION COMPLETE")
        println("="^70)
        println("Collected $(length(snapshots)) snapshots")
        println("\nSnapshot Timeline:")
        for (i, snap) in enumerate(snapshots)
            if i <= 5 || i > length(snapshots) - 5
                @printf("  %2d: t = %.4f, step = %d\n", i, snap.t, snap.step)
            elseif i == 6
                println("  ...")
            end
        end
        println("="^70)
        
        # Launch interactive viewer
        println("\nLaunching Interactive Time-Series Viewer...")
        println("\nViewer Controls:")
        println("  • Time slider: Step through snapshots")
        println("  • Play/Pause/Reset: Animate the time evolution")
        println("  • Quantity buttons: Switch between Density, U, V, W velocities")
        println("  • Isosurface sliders: Adjust visualization levels")
        println("  • Mouse: Rotate (drag), Zoom (scroll)")
        println("="^70)
        
        try
            interactive_3d_timeseries(snapshots, grid, params_with_ic)
        catch e
            @warn "Viewer failed" exception=(e, catch_backtrace())
            println("Trying fallback single-frame viewer...")
            try
                interactive_3d_volume(snapshots[end].M, grid, params_with_ic)
            catch e2
                @warn "Fallback viewer also failed"
            end
        end
    end
else
    # Without snapshots
    if rank == 0
        println("\nRunning without snapshots (snapshot_interval=0)...")
    end
    
    M_final, final_time, time_steps, grid = simulation_runner(params_with_ic)
    
    if rank == 0
        println("\n" * "="^70)
        println("SIMULATION COMPLETE")
        println("="^70)
        println("Final time: $final_time")
        println("Time steps: $time_steps")
        println("="^70)
    end
end

MPI.Finalize()

if rank == 0
    println("\nDone!")
end

