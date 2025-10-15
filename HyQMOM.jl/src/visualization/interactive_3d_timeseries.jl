"""
Interactive 3D Time-Series Visualization

This viewer allows stepping through simulation snapshots over time.
"""

import GLMakie
using Printf

"""
    interactive_3d_timeseries(snapshots, grid, params; kwargs...)

Launch an interactive 3D viewer with time-stepping capability.

# Arguments
- `snapshots`: Vector of (M=M_array, t=time, step=step_num) named tuples
- `grid`: Grid structure with xm, ym, zm
- `params`: Simulation parameters

# Keyword Arguments
- `n_streamlines::Int=8`: Number of streamline seeds
- `vector_step::Int=4`: Subsampling for vector field
- `iso_levels::Vector{Float64}=[0.3, 0.5, 0.7]`: Isosurface levels
"""
function interactive_3d_timeseries(snapshots, grid, params;
                                   n_streamlines=8,
                                   vector_step=4,
                                   streamline_length=50,
                                   iso_levels=[0.3, 0.5, 0.7])
    
    println("\n" * "="^70)
    println("3D TIME-SERIES VIEWER")
    println("="^70)
    println("Features:")
    println("  • Step through simulation snapshots over time")
    println("  • TRUE 3D isosurface contours")
    println("  • Velocity isosurfaces: Blue=positive, Red=negative")
    println("  • Time slider to navigate through evolution")
    println("="^70)
    println("Loaded $(length(snapshots)) snapshots")
    println("  Time range: $(snapshots[1].t) to $(snapshots[end].t)")
    println("="^70)
    
    # Extract grid
    Nx = params.Nx
    Ny = params.Ny
    Nz = params.Nz
    xm = collect(grid.xm)
    ym = collect(grid.ym)
    zm = collect(grid.zm)
    
    # Create figure
    fig = GLMakie.Figure(size=(1800, 1000))
    
    # Main 3D axis
    ax = GLMakie.Axis3(fig[1:3, 1:2], 
                       xlabel="x", ylabel="y", zlabel="z",
                       title="3D Time-Series - Crossing Jets",
                       aspect=:data,
                       azimuth=0.3π,
                       elevation=π/8)
    
    # Control panel
    controls = fig[1:3, 3] = GLMakie.GridLayout()
    
    # Current snapshot index (observable)
    current_snapshot_idx = GLMakie.Observable(1)
    
    # Current quantity
    current_quantity = GLMakie.Observable("Density")
    
    # Quantity buttons
    controls[1, 1] = GLMakie.Label(fig, "Quantity:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    btn_density = GLMakie.Button(fig, label="Density")
    btn_speed = GLMakie.Button(fig, label="Speed")
    btn_u = GLMakie.Button(fig, label="U")
    btn_v = GLMakie.Button(fig, label="V")
    btn_w = GLMakie.Button(fig, label="W")
    btn_pressure = GLMakie.Button(fig, label="Pressure")
    btn_temp = GLMakie.Button(fig, label="Temperature")
    
    controls[2, 1] = GLMakie.hgrid!(btn_density, btn_speed; tellwidth=false)
    controls[3, 1] = GLMakie.hgrid!(btn_u, btn_v, btn_w; tellwidth=false)
    controls[4, 1] = GLMakie.hgrid!(btn_pressure, btn_temp; tellwidth=false)
    
    GLMakie.on(btn_density.clicks) do _
        current_quantity[] = "Density"
    end
    GLMakie.on(btn_speed.clicks) do _
        current_quantity[] = "Speed"
    end
    GLMakie.on(btn_u.clicks) do _
        current_quantity[] = "U velocity"
    end
    GLMakie.on(btn_v.clicks) do _
        current_quantity[] = "V velocity"
    end
    GLMakie.on(btn_w.clicks) do _
        current_quantity[] = "W velocity"
    end
    GLMakie.on(btn_pressure.clicks) do _
        current_quantity[] = "Pressure"
    end
    GLMakie.on(btn_temp.clicks) do _
        current_quantity[] = "Temperature"
    end
    
    # Time controls
    controls[5, 1] = GLMakie.Label(fig, "Time Control:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    # Time slider
    time_slider = GLMakie.Slider(fig, range=1:length(snapshots), startvalue=1)
    controls[6, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Snapshot:", fontsize=12, halign=:left),
        time_slider,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%d/%d", $(time_slider.value), length(snapshots))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Time info label
    time_info = GLMakie.@lift begin
        idx = $(time_slider.value)
        t = snapshots[idx].t
        step = snapshots[idx].step
        @sprintf("t=%.4f, step=%d", t, step)
    end
    controls[7, 1] = GLMakie.Label(fig, time_info, fontsize=12, halign=:left, tellwidth=false)
    
    # Play/Pause buttons (using ASCII to avoid font issues)
    btn_play = GLMakie.Button(fig, label="> Play")
    btn_pause = GLMakie.Button(fig, label="|| Pause")
    btn_reset = GLMakie.Button(fig, label="<< Reset")
    controls[8, 1] = GLMakie.hgrid!(btn_play, btn_pause, btn_reset; tellwidth=false)
    
    is_playing = GLMakie.Observable(false)
    
    GLMakie.on(btn_play.clicks) do _
        is_playing[] = true
    end
    GLMakie.on(btn_pause.clicks) do _
        is_playing[] = false
    end
    GLMakie.on(btn_reset.clicks) do _
        is_playing[] = false
        time_slider.value[] = 1
    end
    
    # Visualization controls
    controls[9, 1] = GLMakie.Label(fig, "Visualization:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    toggle_isosurface = GLMakie.Toggle(fig, active=true)
    controls[10, 1] = GLMakie.hgrid!(
        toggle_isosurface,
        GLMakie.Label(fig, "Isosurfaces", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    toggle_streamlines = GLMakie.Toggle(fig, active=false)
    controls[11, 1] = GLMakie.hgrid!(
        toggle_streamlines,
        GLMakie.Label(fig, "Streamlines", fontsize=12, halign=:left);
        tellwidth=false
    )
    
    # Info about isosurfaces
    controls[12, 1] = GLMakie.Label(fig, 
        "For velocities: Blue=positive, Red=negative", 
        fontsize=10, halign=:left, tellwidth=false)
    
    # Isosurface level sliders
    slider_iso1 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[1])
    controls[13, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso level 1:", fontsize=12, halign=:left),
        slider_iso1,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.2f", $(slider_iso1.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    slider_iso2 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[2])
    controls[14, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso level 2:", fontsize=12, halign=:left),
        slider_iso2,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.2f", $(slider_iso2.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    slider_iso3 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[3])
    controls[15, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso level 3:", fontsize=12, halign=:left),
        slider_iso3,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.2f", $(slider_iso3.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Opacity control
    slider_alpha = GLMakie.Slider(fig, range=0.1:0.1:1.0, startvalue=0.6)
    controls[16, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Transparency:", fontsize=12, halign=:left),
        slider_alpha,
        GLMakie.Label(fig, GLMakie.@lift(@sprintf("%.1f", $(slider_alpha.value))), 
                     fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Compute all quantities for all snapshots (observables will select which to show)
    function compute_quantities(M_snapshot)
        rho = M_snapshot[:, :, :, 1]
        U = M_snapshot[:, :, :, 2] ./ rho
        V = M_snapshot[:, :, :, 6] ./ rho
        W = M_snapshot[:, :, :, 16] ./ rho
        
        speed = sqrt.(U.^2 .+ V.^2 .+ W.^2)
        C200 = M_snapshot[:, :, :, 3] ./ rho .- U.^2
        C020 = M_snapshot[:, :, :, 7] ./ rho .- V.^2
        C002 = M_snapshot[:, :, :, 17] ./ rho .- W.^2
        temperature = (C200 .+ C020 .+ C002) ./ 3.0
        pressure = rho .* temperature
        
        return (rho=rho, U=U, V=V, W=W, speed=speed, pressure=pressure, temperature=temperature)
    end
    
    # Observable for current data
    current_data_obs = GLMakie.@lift begin
        idx = $(time_slider.value)
        q = $(current_quantity)
        quants = compute_quantities(snapshots[idx].M)
        
        if q == "Density"
            quants.rho
        elseif q == "Speed"
            quants.speed
        elseif q == "U velocity"
            quants.U
        elseif q == "V velocity"
            quants.V
        elseif q == "W velocity"
            quants.W
        elseif q == "Pressure"
            quants.pressure
        else # Temperature
            quants.temperature
        end
    end
    
    # Storage for plots
    iso_plots = []
    streamline_plots = []
    
    # Function to create isosurfaces
    function create_isosurfaces!()
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
            # For velocities, check if we have significant values
            pos_level1 = slider_iso1.value[] * data_absmax
            pos_level2 = slider_iso2.value[] * data_absmax
            pos_level3 = slider_iso3.value[] * data_absmax
            
            neg_level1 = -slider_iso1.value[] * data_absmax
            neg_level2 = -slider_iso2.value[] * data_absmax
            neg_level3 = -slider_iso3.value[] * data_absmax
            
            levels = [pos_level1, pos_level2, pos_level3, neg_level1, neg_level2, neg_level3]
            colors = [:blue, :cyan, :lightblue, :red, :orange, :pink]
            alphas = [0.6, 0.5, 0.4, 0.6, 0.5, 0.4] .* slider_alpha.value[]
        else
            level1 = data_min + slider_iso1.value[] * data_range
            level2 = data_min + slider_iso2.value[] * data_range
            level3 = data_min + slider_iso3.value[] * data_range
            
            levels = [level1, level2, level3]
            colors = [:blue, :green, :red]
            alphas = [0.4, 0.5, 0.6] .* slider_alpha.value[]
        end
        
        # Filter out invalid levels (too close to zero or outside data range)
        for (level, color, alpha) in zip(levels, colors, alphas)
            # Skip if level is essentially zero or outside reasonable range
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
    
    # Function to create streamlines
    function create_streamlines!()
        for plot in streamline_plots
            try
                delete!(ax, plot)
            catch
            end
        end
        empty!(streamline_plots)
        
        if !toggle_streamlines.active[]
            return
        end
        
        idx = time_slider.value[]
        quants = compute_quantities(snapshots[idx].M)
        
        # Subsample for streamlines
        step_x = max(1, div(Nx, n_streamlines))
        step_y = max(1, div(Ny, n_streamlines))
        seed_points = []
        for i in 1:step_x:Nx
            for j in 1:step_y:Ny
                for k in 1:max(1, div(Nz,2)):Nz
                    push!(seed_points, GLMakie.Point3f(xm[i], ym[j], zm[k]))
                end
            end
        end
        
        U_interp = quants.U
        V_interp = quants.V
        W_interp = quants.W
        
        try
            for seed in seed_points
                p = GLMakie.streamplot!(ax, 
                    (x, y, z) -> GLMakie.Point3f(
                        # Simple interpolation - in production use proper interpolation
                        0.1, 0.1, 0.0
                    ),
                    xm[1]..xm[end], ym[1]..ym[end], zm[1]..zm[end],
                    colormap=:viridis,
                    linewidth=1.5,
                    alpha=0.7)
                push!(streamline_plots, p)
            end
        catch e
            @warn "Streamlines failed" exception=e
        end
    end
    
    # Initial plots
    create_isosurfaces!()
    create_streamlines!()
    
    # Update plots when time slider changes
    GLMakie.on(time_slider.value) do val
        create_isosurfaces!()
        create_streamlines!()
    end
    
    # Update plots when quantity changes
    GLMakie.on(current_quantity) do q
        create_isosurfaces!()
    end
    
    # Update plots when sliders change
    for slider in [slider_iso1, slider_iso2, slider_iso3, slider_alpha]
        GLMakie.on(slider.value) do val
            create_isosurfaces!()
        end
    end
    
    # Update plots when toggles change
    GLMakie.on(toggle_isosurface.active) do val
        create_isosurfaces!()
    end
    GLMakie.on(toggle_streamlines.active) do val
        create_streamlines!()
    end
    
    # Animation loop for playback
    GLMakie.on(is_playing) do playing
        if playing
            @async begin
                while is_playing[] && time_slider.value[] < length(snapshots)
                    sleep(0.1)  # Adjust playback speed
                    if is_playing[]
                        time_slider.value[] = time_slider.value[] + 1
                    end
                end
                is_playing[] = false
            end
        end
    end
    
    # Display the figure
    display(fig)
    
    println("\n" * "="^70)
    println("TIME-SERIES VIEWER READY!")
    println("="^70)
    println("Controls:")
    println("  • Use time slider to step through snapshots")
    println("  • Click ▶ Play to animate")
    println("  • Click buttons to switch quantities")
    println("  • Adjust iso level sliders for different contour levels")
    println("  • Mouse: drag to rotate, scroll to zoom")
    println("\nPress Enter in terminal to close.")
    println("="^70)
    
    readline()
    
    return fig
end

