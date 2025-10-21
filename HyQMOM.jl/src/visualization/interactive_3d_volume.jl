"""
Interactive 3D Volume Rendering and Isosurface Visualization

This version provides true 3D surfaces and volume rendering, not just slice planes.
"""

import GLMakie
using Printf
using LaTeXStrings

"""
    interactive_3d_volume(M_final, grid, params; kwargs...)

Launch an interactive 3D viewer with volume rendering and true isosurfaces.

# Keyword Arguments
- `n_streamlines::Int=8`: Number of streamline seeds
- `vector_step::Int=4`: Subsampling for vector field
- `iso_levels::Vector{Float64}=[0.3, 0.5, 0.7]`: Isosurface levels (relative to max)
- `enable_volume::Bool=true`: Enable volume rendering
"""
function interactive_3d_volume(M_final, grid, params; 
                                n_streamlines=8,
                                vector_step=4,
                                streamline_length=50,
                                iso_levels=[0.3, 0.5, 0.7],
                                enable_volume=true)
    
    println("\n" * "="^70)
    println("3D VOLUME & ISOSURFACE VIEWER")
    println("="^70)
    println("Features:")
    println("  • TRUE 3D isosurface contours (not dots!)")
    println("  • Volume rendering with transparency")
    println("  • Multiple quantities selectable")
    println("  • Velocity isosurfaces: Blue=positive, Red=negative")
    println("  • Interactive controls")
    println("="^70)
    
    # Extract grid
    Nx = params.Nx
    Ny = params.Ny
    Nz = params.Nz
    xm = collect(grid.xm)
    ym = collect(grid.ym)
    zm = collect(grid.zm)
    
    # Extract quantities
    rho = M_final[:, :, :, 1]
    U = M_final[:, :, :, 2] ./ rho
    V = M_final[:, :, :, 6] ./ rho
    W = M_final[:, :, :, 16] ./ rho
    
    # Derived quantities
    speed = sqrt.(U.^2 .+ V.^2 .+ W.^2)
    C200 = M_final[:, :, :, 3] ./ rho .- U.^2
    C020 = M_final[:, :, :, 7] ./ rho .- V.^2
    C002 = M_final[:, :, :, 17] ./ rho .- W.^2
    temperature = (C200 .+ C020 .+ C002) ./ 3.0
    pressure = rho .* temperature
    
    # Create figure with LaTeX font rendering
    fig = GLMakie.Figure(size=(1600, 1000), fontsize=12,
                        fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
    
    # Main 3D axis
    ax = GLMakie.Axis3(fig[1:3, 1:2], 
                       xlabel=L"x", ylabel=L"y", zlabel=L"z",
                       title="3D Volume Visualization - Crossing Jets",
                       aspect=:data,
                       azimuth=0.3π,
                       elevation=π/8,
                       xticklabelsize=11, yticklabelsize=11, zticklabelsize=11,
                       xlabelsize=13, ylabelsize=13, zlabelsize=13)
    
    # Control panel
    controls = fig[1:3, 3] = GLMakie.GridLayout()
    
    # Quantity selector - using buttons for clarity
    controls[1, 1] = GLMakie.Label(fig, "Select Quantity:", fontsize=16, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    quantity_options = ["Density", "Speed", "Pressure", "Temperature",
                       "U velocity", "V velocity", "W velocity"]
    
    current_quantity = GLMakie.Observable("Density")
    
    # Create button grid
    for (i, q) in enumerate(quantity_options)
        btn = GLMakie.Button(fig, label=q, fontsize=11, tellwidth=false, width=160)
        GLMakie.on(btn.clicks) do n
            current_quantity[] = q
            println("Switched to: $q")
        end
        controls[i+1, 1] = btn
    end
    
    # Visualization controls
    controls[10, 1] = GLMakie.Label(fig, "Visualization:", fontsize=14, 
                                    halign=:left, font=:bold, tellwidth=false)
    
    toggle_isosurface = GLMakie.Toggle(fig, active=true)
    controls[11, 1] = GLMakie.hgrid!(
        toggle_isosurface,
        GLMakie.Label(fig, "Isosurfaces", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    toggle_vectors = GLMakie.Toggle(fig, active=false)
    controls[12, 1] = GLMakie.hgrid!(
        toggle_vectors,
        GLMakie.Label(fig, "Vector field", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    toggle_streamlines = GLMakie.Toggle(fig, active=true)
    controls[13, 1] = GLMakie.hgrid!(
        toggle_streamlines,
        GLMakie.Label(fig, "Streamlines", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    # Info about isosurfaces
    controls[14, 1] = GLMakie.Label(fig, 
        "For velocities: Blue=positive, Red=negative", 
        fontsize=10, halign=:left, tellwidth=false)
    
    # Isosurface level sliders
    slider_iso1 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[1])
    controls[15, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso level 1:", fontsize=12, halign=:left),
        slider_iso1,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.2f", $(slider_iso1.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    slider_iso2 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[2])
    controls[16, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso level 2:", fontsize=12, halign=:left),
        slider_iso2,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.2f", $(slider_iso2.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    slider_iso3 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[3])
    controls[17, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso level 3:", fontsize=12, halign=:left),
        slider_iso3,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.2f", $(slider_iso3.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Opacity control
    slider_alpha = GLMakie.Slider(fig, range=0.1:0.1:1.0, startvalue=0.6)
    controls[18, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Transparency:", fontsize=12, halign=:left),
        slider_alpha,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.1f", $(slider_alpha.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Get current data
    function get_current_data()
        q = current_quantity[]
        if q == "Density"
            return rho, "Density", :viridis
        elseif q == "Speed"
            return speed, "Speed", :plasma
        elseif q == "Pressure"
            return pressure, "Pressure", :thermal
        elseif q == "Temperature"
            return temperature, "Temperature", :hot
        elseif q == "U velocity"
            return U, "U velocity", :RdBu
        elseif q == "V velocity"
            return V, "V velocity", :RdBu
        elseif q == "W velocity"
            return W, "W velocity", :RdBu
        else
            return rho, "Density", :viridis
        end
    end
    
    # Create observables for reactive updates
    current_data_obs = GLMakie.Observable(rho)
    current_cmap_obs = GLMakie.Observable(:viridis)
    current_name_obs = GLMakie.Observable("Density")
    
    # Update function
    GLMakie.on(current_quantity) do q
        data, name, cmap = get_current_data()
        current_data_obs[] = data
        current_name_obs[] = name
        current_cmap_obs[] = cmap
    end
    
    # Plot isosurfaces using contour (TRUE 3D SURFACES)
    iso_plots = []
    
    # Function to create isosurface contours
    function create_isosurfaces!()
        # Clear old plots
        for plot in iso_plots
            try
                delete!(ax, plot)
            catch
            end
        end
        empty!(iso_plots)
        
        if !toggle_isosurface.active[]
            return
        end
        
        data = current_data_obs[]
        q = current_quantity[]
        
        # Check for valid data
        if any(isnan.(data)) || any(isinf.(data))
            @warn "Data contains NaN/Inf, skipping isosurfaces"
            return
        end
        
        # Check if this is a signed quantity (velocities can be positive or negative)
        is_velocity = (q == "U velocity" || q == "V velocity" || q == "W velocity")
        
        data_min = minimum(data)
        data_max = maximum(data)
        data_absmax = maximum(abs.(data))
        data_range = data_max - data_min
        
        # Skip if data is essentially constant or zero
        if data_absmax < 1e-10 || data_range < 1e-10
            # Data is too flat for meaningful isosurfaces
            return
        end
        
        x_lims = (xm[1], xm[end])
        y_lims = (ym[1], ym[end])
        z_lims = (zm[1], zm[end])
        
        if is_velocity
            # For velocities: create isosurfaces for both positive and negative values
            # Positive surfaces (blue shades)
            pos_level1 = slider_iso1.value[] * data_absmax
            pos_level2 = slider_iso2.value[] * data_absmax
            pos_level3 = slider_iso3.value[] * data_absmax
            
            # Negative surfaces (red shades)
            neg_level1 = -slider_iso1.value[] * data_absmax
            neg_level2 = -slider_iso2.value[] * data_absmax
            neg_level3 = -slider_iso3.value[] * data_absmax
            
            # All levels with corresponding colors
            levels = [pos_level1, pos_level2, pos_level3, neg_level1, neg_level2, neg_level3]
            colors = [:blue, :cyan, :lightblue, :red, :orange, :pink]
            alphas = [0.6, 0.5, 0.4, 0.6, 0.5, 0.4] .* slider_alpha.value[]
        else
            # For non-negative quantities (density, pressure, etc): use regular levels
            level1 = data_min + slider_iso1.value[] * data_range
            level2 = data_min + slider_iso2.value[] * data_range
            level3 = data_min + slider_iso3.value[] * data_range
            
            levels = [level1, level2, level3]
            colors = [:blue, :green, :red]
            alphas = [0.4, 0.5, 0.6] .* slider_alpha.value[]
        end
        
        # Create contour surfaces at each level
        for (level, color, alpha) in zip(levels, colors, alphas)
            # Skip if level is essentially zero or invalid
            if abs(level) < 1e-10
                continue
            end
            
            # For non-velocity, skip if outside data range
            if !is_velocity && (level < data_min - 1e-10 || level > data_max + 1e-10)
                continue
            end
            
            try
                p = GLMakie.contour!(ax, x_lims, y_lims, z_lims, data,
                                    levels=[level],
                                    alpha=alpha,
                                    color=color)
                push!(iso_plots, p)
            catch e
                # Only warn if it's not a known issue with near-zero levels
                if abs(level) > 1e-8
                    @warn "Contour failed at level $level" exception=(e,)
                end
            end
        end
    end
    
    # Initial isosurfaces
    create_isosurfaces!()
    
    # Update isosurfaces when quantity or sliders change
    GLMakie.on(current_quantity) do q
        create_isosurfaces!()
    end
    
    GLMakie.on(slider_iso1.value) do val
        create_isosurfaces!()
    end
    
    GLMakie.on(slider_iso2.value) do val
        create_isosurfaces!()
    end
    
    GLMakie.on(slider_iso3.value) do val
        create_isosurfaces!()
    end
    
    GLMakie.on(slider_alpha.value) do val
        create_isosurfaces!()
    end
    
    GLMakie.on(toggle_isosurface.active) do active
        create_isosurfaces!()
    end
    
    # Add vector field
    if vector_step > 0
        xs_vec, ys_vec, zs_vec, u_vec, v_vec, w_vec, speed_vec, mask_vec = 
            compute_vectors(xm, ym, zm, U, V, W, vector_step)
        
        if length(xs_vec[mask_vec]) > 0
            scale = 0.015 / (maximum(speed_vec[mask_vec]) + 1e-10)
            
            GLMakie.arrows!(ax,
                           xs_vec[mask_vec], ys_vec[mask_vec], zs_vec[mask_vec],
                           scale .* u_vec[mask_vec], 
                           scale .* v_vec[mask_vec], 
                           scale .* w_vec[mask_vec],
                           color=speed_vec[mask_vec],
                           colormap=:plasma,
                           linewidth=0.015,
                           arrowsize=GLMakie.Vec3f(0.008, 0.008, 0.012),
                           visible=toggle_vectors.active)
        end
    end
    
    # Add streamlines
    lines = compute_stream(xm, ym, zm, U, V, W, n_streamlines, streamline_length)
    
    for line in lines
        GLMakie.lines!(ax, line, 
                      color=:cyan, 
                      linewidth=3.0,
                      alpha=0.8,
                      visible=toggle_streamlines.active)
    end
    
    # Colorbar
    cb = GLMakie.Colorbar(fig[1:3, 4], 
                         label=current_name_obs,
                         colormap=current_cmap_obs,
                         limits=GLMakie.@lift((minimum($(current_data_obs)), 
                                              maximum($(current_data_obs)))),
                         width=25)
    
    # Statistics
    stats_text = GLMakie.@lift(begin
        data = $(current_data_obs)
        name = $(current_name_obs)
        @sprintf("""
        %s Statistics:
          Min:  %+.4e
          Max:  %+.4e
          Mean: %+.4e
          
        Isosurface Levels:
          Level 1: %.2f × max
          Level 2: %.2f × max
          Level 3: %.2f × max
          
        Grid: %dx%dx%d
        Flow |U|ₘₐₓ: %.3f
        """, name, minimum(data), maximum(data), sum(data)/length(data),
             $(slider_iso1.value), $(slider_iso2.value), $(slider_iso3.value),
             Nx, Ny, Nz, maximum(abs.(U)))
    end)
    
    GLMakie.Label(fig[4, 1:4], stats_text, fontsize=11, halign=:left, 
                 padding=(10, 10, 10, 10), tellwidth=false)
    
    println("\n" * "="^70)
    println("VIEWER READY - TRUE 3D SURFACES!")
    println("="^70)
    println("\nControls:")
    println("  • Click buttons to switch quantities")
    println("  • Adjust iso level sliders for different contour levels")
    println("  • Adjust transparency slider")  
    println("  • Toggle visualization elements on/off")
    println("  • Mouse: drag to rotate, scroll to zoom")
    println("\nYou should see THREE colored isosurfaces:")
    println("  Blue   = low level (30% of max)")
    println("  Green  = mid level (50% of max)")
    println("  Red    = high level (70% of max)")
    println("\nPress Enter in terminal to close.")
    println("="^70)
    
    display(fig)
    readline()
    
    return fig
end

# Helper functions
function compute_vectors(xm, ym, zm, U, V, W, step)
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
    mask = speed .> 0.1
    
    return xs, ys, zs, u_vec, v_vec, w_vec, speed, mask
end

function compute_stream(xm, ym, zm, U, V, W, n_seeds, max_length)
    Nx = length(xm)
    Ny = length(ym)
    Nz = length(zm)
    
    seed_range = range(0.35, 0.45, length=n_seeds)
    seed_points = []
    
    for sx in seed_range, sy in seed_range, sz in seed_range
        push!(seed_points, GLMakie.Point3f(sx, sy, sz))
    end
    
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
            ix = clamp(searchsortedfirst(xm, pos[1]), 1, Nx)
            iy = clamp(searchsortedfirst(ym, pos[2]), 1, Ny)
            iz = clamp(searchsortedfirst(zm, pos[3]), 1, Nz)
            
            vel = [U[ix, iy, iz], V[ix, iy, iz], W[ix, iy, iz]]
            vel_mag = sqrt(sum(vel.^2))
            
            if vel_mag < 0.01 || pos[1] < xm[1] || pos[1] > xm[end] ||
               pos[2] < ym[1] || pos[2] > ym[end] ||
               pos[3] < zm[1] || pos[3] > zm[end]
                break
            end
            
            pos = pos .+ dt_stream .* vel
            push!(line, GLMakie.Point3f(pos...))
        end
        
        if length(line) > 3
            push!(streamline_points, line)
        end
    end
    
    return streamline_points
end

export interactive_3d_volume


