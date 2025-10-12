"""
Simple Interactive 3D Viewer

Launch a basic interactive 3D visualization of the crossing jets.
Uses GLMakie for interactive rotation, zoom, and exploration.

Usage:
    julia --project=. examples/simple_3d_interactive.jl
"""

using HyQMOM
using MPI
using GLMakie
using Printf

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Run simulation (small for quick testing)
params = (
    Np = 40,
    Nz = 20,
    tmax = 0.05,  # Short simulation
    Kn = 1.0,
    Ma = 1.0,
    flag2D = 0,
    CFL = 0.7,
    dx = 1.0/40,
    dy = 1.0/40,
    dz = 1.0/20,
    Nmom = 35,
    nnmax = 100,
    dtmax = 1e-2,
    rhol = 1.0,
    rhor = 0.01,
    T = 1.0,
    r110 = 0.0,
    r101 = 0.0,
    r011 = 0.0,
    symmetry_check_interval = 100,
    homogeneous_z = false,  # TRUE 3D
    debug_output = false,
    enable_memory_tracking = false
)

if rank == 0
    println("Running 3D simulation...")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

# Interactive visualization (only on rank 0)
if rank == 0 && M_final !== nothing
    println("\n" * "="^70)
    println("Launching Interactive 3D Viewer")
    println("="^70)
    println("Controls:")
    println("  - Left mouse drag: Rotate")
    println("  - Right mouse drag: Zoom")
    println("  - Scroll wheel: Zoom")
    println("  - Click 'Density', 'U', 'V', or 'W' to change quantity")
    println("="^70)
    
    # Extract data
    Np = params.Np
    Nz = params.Nz
    xm = grid.xm
    ym = grid.ym
    zm = grid.zm
    
    # Create figure
    fig = GLMakie.Figure(size=(1400, 900))
    
    # Main 3D axis
    ax = GLMakie.Axis3(fig[1, 1:3], 
                       xlabel="x", ylabel="y", zlabel="z",
                       title="3D Crossing Jets - Density",
                       aspect=:data)
    
    # Extract density (avoid naming conflict with GLMakie.density function)
    rho_3d = M_final[:, :, :, 1]
    
    # Create 3D grid of points
    x_pts = Float64[]
    y_pts = Float64[]
    z_pts = Float64[]
    rho_vals = Float64[]
    
    # Sample points where density > threshold
    threshold = 0.1
    for k in 1:Nz, j in 1:Np, i in 1:Np
        if rho_3d[i, j, k] > threshold
            push!(x_pts, xm[i])
            push!(y_pts, ym[j])
            push!(z_pts, zm[k])
            push!(rho_vals, rho_3d[i, j, k])
        end
    end
    
    # Plot as 3D scatter
    GLMakie.scatter!(ax, x_pts, y_pts, z_pts,
                    color=rho_vals,
                    colormap=:viridis,
                    markersize=8,
                    alpha=0.6)
    
    # Add slice planes
    i_mid = div(Np, 2)
    j_mid = div(Np, 2)
    k_mid = div(Nz, 2)
    
    # XY plane
    x_grid = repeat(xm, 1, Np)
    y_grid = repeat(ym', Np, 1)
    z_const = fill(zm[k_mid], Np, Np)
    rho_xy = rho_3d[:, :, k_mid]
    
    GLMakie.surface!(ax, x_grid, y_grid, z_const,
                    color=rho_xy,
                    colormap=:viridis,
                    alpha=0.5)
    
    # Add colorbar
    GLMakie.Colorbar(fig[1, 4], label="Density", colormap=:viridis,
                    limits=(minimum(rho_3d), maximum(rho_3d)))
    
    # Add statistics
    stats_text = @sprintf("""
    Statistics:
      Min:  %.4f
      Max:  %.4f  
      Mean: %.4f
    """, minimum(rho_3d), maximum(rho_3d), sum(rho_3d)/length(rho_3d))
    
    GLMakie.Label(fig[2, 1:4], stats_text, fontsize=14, halign=:left)
    
    println("\nDisplaying interactive viewer...")
    println("Press Enter in terminal to close and exit.")
    
    # Display the figure
    display(fig)
    
    # Keep window open - wait for user input
    readline()
end

MPI.Finalize()

