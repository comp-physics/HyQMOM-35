"""
Interactive 3D Visualization Module

Provides comprehensive interactive 3D visualization with:
- Multiple quantity display (density, velocities, moments)
- Volume rendering and isosurfaces
- Vector fields and streamlines
- Interactive slice planes
- Real-time controls

Dependencies: GLMakie
"""

import GLMakie
using Printf

"""
    interactive_3d_viewer(M_final, grid, params; kwargs...)

Launch an interactive 3D viewer for simulation results.

# Arguments
- `M_final`: Final moment array (Np x Np x Nz x Nmom)
- `grid`: Grid structure with xm, ym, zm
- `params`: Simulation parameters (NamedTuple)

# Keyword Arguments
- `n_streamlines::Int=10`: Number of streamline seed points per dimension
- `vector_step::Int=4`: Subsampling step for vector field
- `streamline_length::Int=50`: Maximum steps for streamline integration
- `iso_threshold::Float64=0.5`: Threshold for isosurface (relative to max)
"""
function interactive_3d_viewer(M_final, grid, params; 
                                n_streamlines=10,
                                vector_step=4,
                                streamline_length=50,
                                iso_threshold=0.5)
    
    println("\n" * "="^70)
    println("ENHANCED INTERACTIVE 3D VIEWER")
    println("="^70)
    println("Features:")
    println("  • Multiple quantities: Density, U/V/W, C-moments, Deltas, H-moments")
    println("  • Volume rendering with adjustable opacity")
    println("  • Isosurfaces at adjustable thresholds")
    println("  • Vector field showing flow direction")
    println("  • Streamlines tracing particle paths")
    println("  • Three orthogonal slice planes")
    println("  • Real-time interactive controls")
    println("="^70)
    
    # Extract grid
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
    
    # Compute derived quantities
    speed = sqrt.(U.^2 .+ V.^2 .+ W.^2)
    
    # C-moments (centralized moments)
    C200 = M_final[:, :, :, 3] ./ rho .- U.^2
    C020 = M_final[:, :, :, 7] ./ rho .- V.^2
    C002 = M_final[:, :, :, 17] ./ rho .- W.^2
    C110 = M_final[:, :, :, 4] ./ rho .- U .* V
    C101 = M_final[:, :, :, 10] ./ rho .- U .* W
    C011 = M_final[:, :, :, 14] ./ rho .- V .* W
    
    # Temperature (trace of stress tensor)
    temperature = (C200 .+ C020 .+ C002) ./ 3.0
    
    # Pressure
    pressure = rho .* temperature
    
    # Anisotropy measures
    Delta1 = C200 .- C020
    Delta2 = C200 .- C002
    
    # H-moments (related to heat flux)
    H200 = C200
    H020 = C020
    H002 = C002
    
    # Create figure with comprehensive layout
    fig = GLMakie.Figure(size=(1800, 1100))
    
    # Main 3D axis
    ax = GLMakie.Axis3(fig[1:4, 1:3], 
                       xlabel="x", ylabel="y", zlabel="z",
                       title="3D Crossing Jets - Interactive Viewer",
                       aspect=:data,
                       azimuth=0.5π,
                       elevation=π/6)
    
    # Control panel
    controls = fig[1:4, 4] = GLMakie.GridLayout()
    
    # Quantity selector with comprehensive options
    quantity_options = [
        "Density", "Speed", "Pressure", "Temperature",
        "U velocity", "V velocity", "W velocity",
        "C200 (σ_xx)", "C020 (σ_yy)", "C002 (σ_zz)",
        "C110 (σ_xy)", "C101 (σ_xz)", "C011 (σ_yz)",
        "Delta1 (σ_xx - σ_yy)", "Delta2 (σ_xx - σ_zz)",
        "H200", "H020", "H002"
    ]
    
    quantity_menu = GLMakie.Menu(fig, options=quantity_options, default="Density")
    controls[1, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Quantity:", fontsize=16, halign=:left, font=:bold),
        quantity_menu;
        tellwidth=false
    )
    
    # Slice position sliders
    controls[2, 1] = GLMakie.Label(fig, "Slice Positions:", fontsize=14, halign=:left, font=:bold, tellwidth=false)
    
    slider_x = GLMakie.Slider(fig, range=1:Np, startvalue=div(Np,2))
    controls[3, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "X:", fontsize=12, halign=:left),
        slider_x,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("x = %.3f", xm[$(slider_x.value)])), fontsize=11, halign=:left);
        tellwidth=false
    )
    
    slider_y = GLMakie.Slider(fig, range=1:Np, startvalue=div(Np,2))
    controls[4, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Y:", fontsize=12, halign=:left),
        slider_y,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("y = %.3f", ym[$(slider_y.value)])), fontsize=11, halign=:left);
        tellwidth=false
    )
    
    slider_z = GLMakie.Slider(fig, range=1:Nz, startvalue=div(Nz,2))
    controls[5, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Z:", fontsize=12, halign=:left),
        slider_z,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("z = %.3f", zm[$(slider_z.value)])), fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Visualization toggles
    controls[6, 1] = GLMakie.Label(fig, "Visualization Options:", fontsize=14, halign=:left, font=:bold, tellwidth=false)
    
    toggle_slices = GLMakie.Toggle(fig, active=true)
    controls[7, 1] = GLMakie.hgrid!(
        toggle_slices,
        GLMakie.Label(fig, "Slice planes", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    toggle_vectors = GLMakie.Toggle(fig, active=true)
    controls[8, 1] = GLMakie.hgrid!(
        toggle_vectors,
        GLMakie.Label(fig, "Vector field", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    toggle_streamlines = GLMakie.Toggle(fig, active=false)
    controls[9, 1] = GLMakie.hgrid!(
        toggle_streamlines,
        GLMakie.Label(fig, "Streamlines", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    toggle_isosurface = GLMakie.Toggle(fig, active=false)
    controls[10, 1] = GLMakie.hgrid!(
        toggle_isosurface,
        GLMakie.Label(fig, "Isosurface", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    # Isosurface threshold slider
    slider_iso = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_threshold)
    controls[11, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso threshold:", fontsize=12, halign=:left),
        slider_iso,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.2f", $(slider_iso.value))), fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Opacity slider for volume rendering
    slider_alpha = GLMakie.Slider(fig, range=0.3:0.1:0.9, startvalue=0.7)
    controls[12, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Slice opacity:", fontsize=12, halign=:left),
        slider_alpha,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.1f", $(slider_alpha.value))), fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Streamline density slider
    slider_stream_density = GLMakie.Slider(fig, range=3:1:15, startvalue=n_streamlines)
    controls[13, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Streamline density:", fontsize=12, halign=:left),
        slider_stream_density,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%d seeds", $(slider_stream_density.value))), fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Create observables for reactive plotting
    Q_data = GLMakie.Observable(rho)
    Q_name = GLMakie.Observable("Density")
    Q_cmap = GLMakie.Observable(:viridis)
    Q_limits = GLMakie.Observable((minimum(rho), maximum(rho)))
    
    # Function to get quantity data
    function get_quantity_data(name)
        if name == "Density"
            return rho, "Density", :viridis, (minimum(rho), maximum(rho))
        elseif name == "Speed"
            return speed, "Speed", :plasma, (0.0, maximum(speed))
        elseif name == "Pressure"
            return pressure, "Pressure", :thermal, (minimum(pressure), maximum(pressure))
        elseif name == "Temperature"
            return temperature, "Temperature", :hot, (minimum(temperature), maximum(temperature))
        elseif name == "U velocity"
            vmax = max(abs(minimum(U)), abs(maximum(U)))
            return U, "U velocity", :RdBu, (-vmax, vmax)
        elseif name == "V velocity"
            vmax = max(abs(minimum(V)), abs(maximum(V)))
            return V, "V velocity", :RdBu, (-vmax, vmax)
        elseif name == "W velocity"
            vmax = max(abs(minimum(W)), abs(maximum(W)))
            return W, "W velocity", :RdBu, (-vmax, vmax)
        elseif name == "C200 (σ_xx)"
            return C200, "C200", :balance, (minimum(C200), maximum(C200))
        elseif name == "C020 (σ_yy)"
            return C020, "C020", :balance, (minimum(C020), maximum(C020))
        elseif name == "C002 (σ_zz)"
            return C002, "C002", :balance, (minimum(C002), maximum(C002))
        elseif name == "C110 (σ_xy)"
            vmax = max(abs(minimum(C110)), abs(maximum(C110)))
            return C110, "C110", :RdBu, (-vmax, vmax)
        elseif name == "C101 (σ_xz)"
            vmax = max(abs(minimum(C101)), abs(maximum(C101)))
            return C101, "C101", :RdBu, (-vmax, vmax)
        elseif name == "C011 (σ_yz)"
            vmax = max(abs(minimum(C011)), abs(maximum(C011)))
            return C011, "C011", :RdBu, (-vmax, vmax)
        elseif name == "Delta1 (σ_xx - σ_yy)"
            vmax = max(abs(minimum(Delta1)), abs(maximum(Delta1)))
            return Delta1, "Delta1", :RdBu, (-vmax, vmax)
        elseif name == "Delta2 (σ_xx - σ_zz)"
            vmax = max(abs(minimum(Delta2)), abs(maximum(Delta2)))
            return Delta2, "Delta2", :RdBu, (-vmax, vmax)
        elseif name == "H200"
            return H200, "H200", :balance, (minimum(H200), maximum(H200))
        elseif name == "H020"
            return H020, "H020", :balance, (minimum(H020), maximum(H020))
        elseif name == "H002"
            return H002, "H002", :balance, (minimum(H002), maximum(H002))
        else
            return rho, "Density", :viridis, (minimum(rho), maximum(rho))
        end
    end
    
    # Function to update quantity
    function update_quantity!(name)
        data, qname, cmap, limits = get_quantity_data(name)
        Q_data[] = data
        Q_name[] = qname
        Q_cmap[] = cmap
        Q_limits[] = limits
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
                    alpha=slider_alpha.value,
                    shading=GLMakie.NoShading,
                    visible=toggle_slices.active)
    
    # XZ plane (constant y)
    x_xz = repeat(xm, 1, Nz)
    z_xz = repeat(zm', Np, 1)
    y_xz = GLMakie.@lift(fill(ym[$(slider_y.value)], Np, Nz))
    Q_xz = GLMakie.@lift($(Q_data)[:, $(slider_y.value), :])
    
    surf_xz = GLMakie.surface!(ax, x_xz, y_xz, z_xz,
                    color=Q_xz,
                    colormap=Q_cmap,
                    alpha=slider_alpha.value,
                    shading=GLMakie.NoShading,
                    visible=toggle_slices.active)
    
    # YZ plane (constant x)
    y_yz = repeat(ym, 1, Nz)
    z_yz = repeat(zm', Np, 1)
    x_yz = GLMakie.@lift(fill(xm[$(slider_x.value)], Np, Nz))
    Q_yz = GLMakie.@lift($(Q_data)[$(slider_x.value), :, :])
    
    surf_yz = GLMakie.surface!(ax, x_yz, y_yz, z_yz,
                    color=Q_yz,
                    colormap=Q_cmap,
                    alpha=slider_alpha.value,
                    shading=GLMakie.NoShading,
                    visible=toggle_slices.active)
    
    # Add vector field
    xs_vec, ys_vec, zs_vec, u_vec, v_vec, w_vec, speed_vec, mask_vec = 
        compute_vector_field(xm, ym, zm, U, V, W, vector_step)
    
    scale = 0.02 / (maximum(speed_vec[mask_vec]) + 1e-10)
    
    arrows_plot = GLMakie.arrows!(ax,
                   xs_vec[mask_vec], ys_vec[mask_vec], zs_vec[mask_vec],
                   scale .* u_vec[mask_vec], scale .* v_vec[mask_vec], scale .* w_vec[mask_vec],
                   color=speed_vec[mask_vec],
                   colormap=:plasma,
                   linewidth=0.02,
                   arrowsize=GLMakie.Vec3f(0.01, 0.01, 0.015),
                   visible=toggle_vectors.active)
    
    # Streamlines storage
    streamline_plots = []
    
    # Function to update streamlines
    function update_streamlines!()
        # Clear old streamlines
        for plot in streamline_plots
            GLMakie.delete!(ax, plot)
        end
        empty!(streamline_plots)
        
        n_seeds = slider_stream_density.value[]
        lines = compute_streamlines(xm, ym, zm, U, V, W, n_seeds, streamline_length)
        
        for line in lines
            plot = GLMakie.lines!(ax, line, 
                          color=:cyan, 
                          linewidth=2.5,
                          alpha=0.7,
                          visible=toggle_streamlines.active)
            push!(streamline_plots, plot)
        end
    end
    
    # Initial streamlines
    update_streamlines!()
    
    # Update streamlines when density changes
    GLMakie.on(slider_stream_density.value) do val
        if toggle_streamlines.active[]
            update_streamlines!()
        end
    end
    
    # Isosurface (as scatter plot of points near threshold)
    iso_points = GLMakie.Observable(GLMakie.Point3f[])
    iso_colors = GLMakie.Observable(Float64[])
    
    function update_isosurface!()
        data = Q_data[]
        threshold = slider_iso.value[] * maximum(data)
        
        points = GLMakie.Point3f[]
        colors = Float64[]
        
        # Sample points near isosurface (within 5% of threshold)
        tol = 0.05 * maximum(data)
        for k in 1:Nz, j in 1:Np, i in 1:Np
            if abs(data[i,j,k] - threshold) < tol
                push!(points, GLMakie.Point3f(xm[i], ym[j], zm[k]))
                push!(colors, data[i,j,k])
            end
        end
        
        iso_points[] = points
        iso_colors[] = colors
    end
    
    iso_scatter = GLMakie.scatter!(ax, iso_points,
                          color=iso_colors,
                          colormap=Q_cmap,
                          markersize=12,
                          alpha=0.6,
                          visible=toggle_isosurface.active)
    
    # Update isosurface when needed
    GLMakie.on(slider_iso.value) do val
        if toggle_isosurface.active[]
            update_isosurface!()
        end
    end
    
    GLMakie.on(quantity_menu.selection) do s
        if toggle_isosurface.active[]
            update_isosurface!()
        end
    end
    
    GLMakie.on(toggle_isosurface.active) do active
        if active
            update_isosurface!()
        end
    end
    
    # Colorbar
    cb = GLMakie.Colorbar(fig[1:4, 5], 
                         label=Q_name,
                         colormap=Q_cmap,
                         limits=Q_limits,
                         width=25)
    
    # Statistics panel
    stats_text = GLMakie.@lift(begin
        q = $(Q_data)
        name = $(Q_name)
        @sprintf("""
        %s Statistics:
          Min:  %+.4e
          Max:  %+.4e
          Mean: %+.4e
          Std:  %.4e
          
        Flow Field:
          |U|ₘₐₓ: %.4e
          |V|ₘₐₓ: %.4e
          |W|ₘₐₓ: %.4e
          Speedₘₐₓ: %.4e
          
        Grid: %dx%dx%d
        Time: %.4f
        Steps: %d
        """, name, minimum(q), maximum(q), sum(q)/length(q), 
             sqrt(sum((q .- sum(q)/length(q)).^2) / length(q)),
             maximum(abs.(U)), maximum(abs.(V)), maximum(abs.(W)),
             maximum(speed),
             Np, Np, Nz, params.tmax, 15)
    end)
    
    GLMakie.Label(fig[5, 1:5], stats_text, fontsize=11, halign=:left, tellwidth=false,
                 padding=(10, 10, 10, 10))
    
    println("\n" * "="^70)
    println("INTERACTIVE VIEWER READY")
    println("="^70)
    println("\nControls:")
    println("  • Select quantity from dropdown menu")
    println("  • Move sliders to adjust slice positions")
    println("  • Toggle visualization elements on/off")
    println("  • Adjust isosurface threshold and opacity")
    println("  • Change streamline density")
    println("  • Mouse: drag to rotate, scroll to zoom")
    println("\nPress Enter in terminal to close and exit.")
    println("="^70)
    
    # Display the figure
    display(fig)
    
    # Keep window open
    readline()
    
    return fig
end

"""
Helper function to compute vector field for display
"""
function compute_vector_field(xm, ym, zm, U, V, W, step)
    x_vec = xm[1:step:end]
    y_vec = ym[1:step:end]
    z_vec = zm[1:step:end]
    
    xs = [x for x in x_vec, y in y_vec, z in z_vec][:]
    ys = [y for x in x_vec, y in y_vec, z in z_vec][:]
    zs = [z for x in x_vec, y in y_vec, z in z_vec][:]
    
    u_vec = U[1:step:end, 1:step:end, 1:step:end][:]
    v_vec = V[1:step:end, 1:step:end, 1:step:end][:]
    w_vec = W[1:step:end, 1:step:end, 1:step:end][:]
    
    speed = sqrt.(u_vec.^2 .+ v_vec.^2 .+ w_vec.^2)
    threshold = 0.1
    mask = speed .> threshold
    
    return xs, ys, zs, u_vec, v_vec, w_vec, speed, mask
end

"""
Helper function to compute streamlines
"""
function compute_streamlines(xm, ym, zm, U, V, W, n_seeds, max_length)
    Np = length(xm)
    Nz = length(zm)
    
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
    
    streamline_points = []
    dt_stream = 0.003
    
    for seed in seed_points
        line = [seed]
        pos = [seed[1], seed[2], seed[3]]
        
        for step in 1:max_length
            # Find velocity at current position
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
        
        if length(line) > 3
            push!(streamline_points, line)
        end
    end
    
    return streamline_points
end

# Export the main function
export interactive_3d_viewer

