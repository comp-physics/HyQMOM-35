"""
Enhanced Interactive 3D Viewer with Flow Visualization

Features:
- Toggle between density, U, V, W velocities
- Vector field showing velocity arrows
- Streamlines for flow visualization
- Adjustable slice planes
- Interactive controls with sliders

Usage:
    julia --project=. examples/enhanced_3d_interactive.jl
"""

using HyQMOM
using MPI
using GLMakie
using Printf

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Run simulation
params = (
    Np = 40,
    Nz = 20,
    tmax = 0.01,
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

# Enhanced interactive visualization
if rank == 0 && M_final !== nothing
    println("\n" * "="^70)
    println("Enhanced Interactive 3D Viewer")
    println("="^70)
    println("Features:")
    println("  - Toggle quantities with menu")
    println("  - Adjust slice positions with sliders")
    println("  - Vector field shows flow direction")
    println("  - Streamlines trace particle paths")
    println("  - Rotate/zoom/pan with mouse")
    println("="^70)
    
    # Extract data
    Np = params.Np
    Nz = params.Nz
    xm = collect(grid.xm)
    ym = collect(grid.ym)
    zm = collect(grid.zm)
    
    # Extract all quantities
    rho = M_final[:, :, :, 1]
    U = M_final[:, :, :, 2] ./ rho
    V = M_final[:, :, :, 6] ./ rho
    W = M_final[:, :, :, 16] ./ rho
    
    # Create figure with controls
    fig = GLMakie.Figure(size=(1600, 1000))
    
    # Main 3D axis
    ax = GLMakie.Axis3(fig[1:3, 1:2], 
                       xlabel="x", ylabel="y", zlabel="z",
                       title="3D Crossing Jets - Interactive Viewer",
                       aspect=:data,
                       azimuth=0.5π,
                       elevation=π/6)
    
    # Control panel
    controls = fig[1:3, 3] = GLMakie.GridLayout()
    
    # Quantity selector
    quantity_menu = GLMakie.Menu(fig, options=["Density", "U velocity", "V velocity", "W velocity"])
    controls[1, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Quantity:", fontsize=16, halign=:left),
        quantity_menu
    )
    
    # Slice position sliders
    slider_x = GLMakie.Slider(fig, range=1:Np, startvalue=div(Np,2))
    slider_y = GLMakie.Slider(fig, range=1:Np, startvalue=div(Np,2))
    slider_z = GLMakie.Slider(fig, range=1:Nz, startvalue=div(Nz,2))
    
    controls[2, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "X slice:", fontsize=14, halign=:left),
        slider_x,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("x = %.3f", xm[$(slider_x.value)])), fontsize=12, halign=:left)
    )
    
    controls[3, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Y slice:", fontsize=14, halign=:left),
        slider_y,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("y = %.3f", ym[$(slider_y.value)])), fontsize=12, halign=:left)
    )
    
    controls[4, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Z slice:", fontsize=14, halign=:left),
        slider_z,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("z = %.3f", zm[$(slider_z.value)])), fontsize=12, halign=:left)
    )
    
    # Toggle for vector field
    toggle_vectors = GLMakie.Toggle(fig, active=true)
    controls[5, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Show vectors:", fontsize=14, halign=:left),
        toggle_vectors
    )
    
    # Toggle for streamlines
    toggle_streamlines = GLMakie.Toggle(fig, active=false)
    controls[6, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Show streamlines:", fontsize=14, halign=:left),
        toggle_streamlines
    )
    
    # Create observables for reactive plotting
    Q_data = GLMakie.Observable(rho)
    Q_name = GLMakie.Observable("Density")
    Q_cmap = GLMakie.Observable(:viridis)
    Q_limits = GLMakie.Observable((minimum(rho), maximum(rho)))
    
    # Function to update quantity when menu changes
    function update_quantity!(name)
        if name == "Density"
            Q_data[] = rho
            Q_name[] = "Density"
            Q_cmap[] = :viridis
            Q_limits[] = (minimum(rho), maximum(rho))
        elseif name == "U velocity"
            Q_data[] = U
            Q_name[] = "U velocity"
            Q_cmap[] = :RdBu
            vmax = max(abs(minimum(U)), abs(maximum(U)))
            Q_limits[] = (-vmax, vmax)
        elseif name == "V velocity"
            Q_data[] = V
            Q_name[] = "V velocity"
            Q_cmap[] = :RdBu
            vmax = max(abs(minimum(V)), abs(maximum(V)))
            Q_limits[] = (-vmax, vmax)
        elseif name == "W velocity"
            Q_data[] = W
            Q_name[] = "W velocity"
            Q_cmap[] = :RdBu
            vmax = max(abs(minimum(W)), abs(maximum(W)))
            Q_limits[] = (-vmax, vmax)
        end
    end
    
    # Connect menu to update function
    GLMakie.on(quantity_menu.selection) do s
        update_quantity!(s)
    end
    
    # XY plane (constant z)
    x_xy = repeat(xm, 1, Np)
    y_xy = repeat(ym', Np, 1)
    z_xy = GLMakie.@lift(fill(zm[$(slider_z.value)], Np, Np))
    Q_xy = GLMakie.@lift($(Q_data)[:, :, $(slider_z.value)])
    
    surf_xy = GLMakie.surface!(ax, x_xy, y_xy, z_xy,
                    color=Q_xy,
                    colormap=Q_cmap,
                    alpha=0.7,
                    shading=GLMakie.NoShading)
    
    # XZ plane (constant y)
    x_xz = repeat(xm, 1, Nz)
    z_xz = repeat(zm', Np, 1)
    y_xz = GLMakie.@lift(fill(ym[$(slider_y.value)], Np, Nz))
    Q_xz = GLMakie.@lift($(Q_data)[:, $(slider_y.value), :])
    
    surf_xz = GLMakie.surface!(ax, x_xz, y_xz, z_xz,
                    color=Q_xz,
                    colormap=Q_cmap,
                    alpha=0.7,
                    shading=GLMakie.NoShading)
    
    # YZ plane (constant x)
    y_yz = repeat(ym, 1, Nz)
    z_yz = repeat(zm', Np, 1)
    x_yz = GLMakie.@lift(fill(xm[$(slider_x.value)], Np, Nz))
    Q_yz = GLMakie.@lift($(Q_data)[$(slider_x.value), :, :])
    
    surf_yz = GLMakie.surface!(ax, x_yz, y_yz, z_yz,
                    color=Q_yz,
                    colormap=Q_cmap,
                    alpha=0.7,
                    shading=GLMakie.NoShading)
    
    # Add vector field (subsample for clarity)
    step = 4
    x_vec = xm[1:step:end]
    y_vec = ym[1:step:end]
    z_vec = zm[1:step:end]
    
    # Create grid for vectors
    xs = [x for x in x_vec, y in y_vec, z in z_vec][:]
    ys = [y for x in x_vec, y in y_vec, z in z_vec][:]
    zs = [z for x in x_vec, y in y_vec, z in z_vec][:]
    
    # Sample velocity at vector grid points
    u_vec = U[1:step:end, 1:step:end, 1:step:end][:]
    v_vec = V[1:step:end, 1:step:end, 1:step:end][:]
    w_vec = W[1:step:end, 1:step:end, 1:step:end][:]
    
    # Scale vectors
    speed = sqrt.(u_vec.^2 .+ v_vec.^2 .+ w_vec.^2)
    scale = 0.02 / (maximum(speed) + 1e-10)
    
    # Only show vectors where speed is significant
    threshold = 0.1
    mask = speed .> threshold
    
    arrows_plot = GLMakie.arrows!(ax,
                   xs[mask], ys[mask], zs[mask],
                   scale .* u_vec[mask], scale .* v_vec[mask], scale .* w_vec[mask],
                   color=speed[mask],
                   colormap=:plasma,
                   linewidth=0.02,
                   arrowsize=GLMakie.Vec3f(0.01, 0.01, 0.015),
                   visible=toggle_vectors.active)
    
    # Add streamlines starting from seed points in the jet regions
    # Seed points in bottom-left and top-right jet regions
    n_seeds = 5
    seed_range = range(0.35, 0.45, length=n_seeds)
    seed_points = []
    
    # Bottom-left jet
    for sx in seed_range, sy in seed_range, sz in seed_range
        push!(seed_points, GLMakie.Point3f(sx, sy, sz))
    end
    
    # Top-right jet  
    for sx in range(0.55, 0.65, length=n_seeds)
        for sy in range(0.55, 0.65, length=n_seeds)
            for sz in seed_range
                push!(seed_points, GLMakie.Point3f(sx, sy, sz))
            end
        end
    end
    
    # Create simple streamlines (particle tracing)
    streamline_points = []
    max_length = 50
    dt_stream = 0.005
    
    for seed in seed_points
        line = [seed]
        pos = [seed[1], seed[2], seed[3]]
        
        for step in 1:max_length
            # Find velocity at current position (simple nearest neighbor)
            ix = clamp(searchsortedfirst(xm, pos[1]), 1, Np)
            iy = clamp(searchsortedfirst(ym, pos[2]), 1, Np)
            iz = clamp(searchsortedfirst(zm, pos[3]), 1, Nz)
            
            vel = [U[ix, iy, iz], V[ix, iy, iz], W[ix, iy, iz]]
            vel_mag = sqrt(sum(vel.^2))
            
            # Stop if velocity too small or out of bounds
            if vel_mag < 0.01 || pos[1] < xm[1] || pos[1] > xm[end] ||
               pos[2] < ym[1] || pos[2] > ym[end] ||
               pos[3] < zm[1] || pos[3] > zm[end]
                break
            end
            
            # Euler step
            pos = pos .+ dt_stream .* vel
            push!(line, GLMakie.Point3f(pos...))
        end
        
        if length(line) > 2
            push!(streamline_points, line)
        end
    end
    
    # Plot streamlines
    for line in streamline_points
        GLMakie.lines!(ax, line, 
                      color=:cyan, 
                      linewidth=2.0,
                      alpha=0.6,
                      visible=toggle_streamlines.active)
    end
    
    # Colorbar
    cb = GLMakie.Colorbar(fig[1:3, 4], 
                         label=Q_name,
                         colormap=Q_cmap,
                         limits=Q_limits)
    
    # Statistics
    stats_text = GLMakie.@lift(begin
        q = $(Q_data)
        name = $(Q_name)
        @sprintf("""
        %s Statistics:
          Min:  %.4f
          Max:  %.4f
          Mean: %.4f
          Std:  %.4f
          
        Flow Stats:
          |U| max: %.4f
          |V| max: %.4f
          |W| max: %.4f
        """, name, minimum(q), maximum(q), sum(q)/length(q), 
             sqrt(sum((q .- sum(q)/length(q)).^2) / length(q)),
             maximum(abs.(U)), maximum(abs.(V)), maximum(abs.(W)))
    end)
    
    GLMakie.Label(fig[4, 1:4], stats_text, fontsize=12, halign=:left, tellwidth=false)
    
    println("\nDisplaying enhanced interactive viewer...")
    println("  - Use the menu to switch between quantities")
    println("  - Adjust sliders to move slice planes")
    println("  - Toggle vector field and streamlines")
    println("  - Drag with mouse to rotate view")
    println("  - Streamlines show particle paths (cyan)")
    println("\nPress Enter in terminal to close and exit.")
    
    # Display the figure
    display(fig)
    
    # Keep window open
    readline()
end

MPI.Finalize()
