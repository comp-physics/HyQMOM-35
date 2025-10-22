"""
Interactive 3D Time-Series Visualization

This viewer allows stepping through simulation snapshots over time.
"""

import GLMakie
using Printf
using LaTeXStrings

# Import moment computation functions
import ..get_standardized_moment

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
    
    # Check if snapshots have standardized moments
    has_std_moments = haskey(snapshots[1], :S)
    
    println("\n" * "="^70)
    println("VIEWER LAYOUT")
    println("="^70)
    println("Has standardized moments (S field)? ", has_std_moments)
    
    # Create figure - default 1800x1000 for both layouts, minimal spacing
    # Use LaTeX rendering for all text (tick labels, axis labels, titles)
    if has_std_moments
        println("Creating 3-column figure (1800×1000):")
        println("  Column 1: Physical space (x,y,z)")
        println("  Column 2: Moment space (S110, S101, S011)")
        println("  Column 3: Controls")
        # Same width as standard, just split among 3 columns, minimal gaps
        fig = GLMakie.Figure(size=(1600, 700), fontsize=12,
                            fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
        # Reduce column gaps to maximize plot space
        GLMakie.colgap!(fig.layout, 5)  # 5 pixels between columns
    else
        println("Creating 2-column figure (1800×1000)")
        # Default size: 1800x1000
        fig = GLMakie.Figure(size=(1600, 700), fontsize=12,
                            fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
        GLMakie.colgap!(fig.layout, 5)  # 5 pixels between columns
    end
    
    # Current snapshot index (observable) - needed for moment space title
    current_snapshot_idx = GLMakie.Observable(1)
    
    # Left: Physical space (isosurfaces)
    ax_physical = GLMakie.Axis3(fig[1, 1], 
                                xlabel=L"x", ylabel=L"y", zlabel=L"z",
                                title="Physical Space - Crossing Jets",
                                aspect=:data,
                                azimuth=0.3π,
                                elevation=π/8,
                                xticklabelsize=11, yticklabelsize=11, zticklabelsize=11,
                                xlabelsize=13, ylabelsize=13, zlabelsize=13)
    
    # Middle: Moment space (if available)
    ax_moment = nothing
    if has_std_moments
        println("Creating moment space axis at position [1, 2]...")
        ax_moment = GLMakie.Axis3(fig[1, 2], 
                                 xlabel=L"S_{110}", ylabel=L"S_{101}", zlabel=L"S_{011}",
                                 title=GLMakie.@lift(latexstring("Moment Space - ", 
                                                                @sprintf("t=%.4f", snapshots[$current_snapshot_idx].t))),
                                 aspect=:data,
                                 azimuth=0.3π,
                                 elevation=π/8,
                                 limits=(-1, 1, -1, 1, -1, 1),
                                 xticklabelsize=11, yticklabelsize=11, zticklabelsize=11,
                                 xlabelsize=13, ylabelsize=13, zlabelsize=13)
        
        println("✓ Moment space axis created successfully!")
        println("  This will show S110, S101, S011 as a 3D scatter plot")
    end
    println("="^70)
    
    # Control panel (far right)
    control_col = has_std_moments ? 3 : 2
    controls = fig[1, control_col] = GLMakie.GridLayout()
    
    # Current quantity
    current_quantity = GLMakie.Observable("Density")
    
    # Quantity buttons (minimal - single row)
    btn_density = GLMakie.Button(fig, label=L"\rho", fontsize=8)
    btn_u = GLMakie.Button(fig, label=L"u", fontsize=8)
    btn_v = GLMakie.Button(fig, label=L"v", fontsize=8)
    btn_w = GLMakie.Button(fig, label=L"w", fontsize=8)
    btn_pressure = GLMakie.Button(fig, label=L"P", fontsize=8)
    
    controls[1, 1] = GLMakie.hgrid!(btn_density, btn_u, btn_v, btn_w, btn_pressure; tellwidth=false)
    
    # Removed: Speed and Temperature buttons to save space
    btn_speed = btn_density  # Map to density
    btn_temp = btn_pressure  # Map to pressure
    
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
    
    # Time slider with snapshot number label
    time_slider = GLMakie.Slider(fig, range=1:length(snapshots), startvalue=1, width=200)
    btn_play = GLMakie.Button(fig, label=">", fontsize=8, width=95)
    btn_pause = GLMakie.Button(fig, label="||", fontsize=8, width=95)
    
    time_label = GLMakie.@lift(@sprintf("Snap %d/%d", $(time_slider.value), length(snapshots)))
    
    controls[2, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, time_label, fontsize=9, halign=:left),
        time_slider;
        tellwidth=false, tellheight=false
    )
    controls[3, 1] = GLMakie.hgrid!(btn_play, btn_pause; tellwidth=false, tellheight=false)
    
    # Add PNG export buttons - separate exports for each plot
    btn_export_physical = GLMakie.Button(fig, label="Phys PNG", fontsize=8, width=95)
    btn_export_moment = GLMakie.Button(fig, label="Mom PNG", fontsize=8, width=95)
    controls[4, 1] = GLMakie.hgrid!(btn_export_physical, btn_export_moment; tellwidth=false, tellheight=false)
    
    is_playing = GLMakie.Observable(false)
    
    GLMakie.on(btn_play.clicks) do _
        is_playing[] = true
    end
    GLMakie.on(btn_pause.clicks) do _
        is_playing[] = false
    end
    
    # PNG export callbacks - separate high-resolution files for each plot
    # Note: We use screen capture to avoid creating new windows that might interfere
    GLMakie.on(btn_export_physical.clicks) do _
        try
            snap_idx = time_slider.value[]
            current_time = snapshots[snap_idx].t
            filename = @sprintf("physical_space_t%.4f_snap%03d.png", current_time, snap_idx)
            
            println("\n" * "="^70)
            println("EXPORTING PHYSICAL SPACE PLOT TO PNG")
            println("="^70)
            println("Filename: $filename")
            println("Snapshot: $snap_idx / $(length(snapshots))")
            println("Time: $current_time")
            println("Capturing current view from screen...")
            
            # Directly save the axis scene to file (captures current view)
            # This avoids creating a new window/figure
            GLMakie.save(filename, ax_physical.blockscene, px_per_unit=2)
            
            println("✓ Physical space plot exported successfully!")
            println("  File size: $(round(filesize(filename)/1024, digits=1)) KB")
            println("="^70)
        catch e
            @warn "Physical space PNG export failed" exception=(e, catch_backtrace())
            println("Error: $e")
        end
    end
    
    GLMakie.on(btn_export_moment.clicks) do _
        try
            snap_idx = time_slider.value[]
            current_time = snapshots[snap_idx].t
            filename = @sprintf("moment_space_t%.4f_snap%03d.png", current_time, snap_idx)
            
            println("\n" * "="^70)
            println("EXPORTING MOMENT SPACE PLOT TO PNG")
            println("="^70)
            println("Filename: $filename")
            println("Snapshot: $snap_idx / $(length(snapshots))")
            println("Time: $current_time")
            
            if !has_std_moments
                println("⚠ No moment space data available - skipping export")
                println("="^70)
                return
            end
            
            println("Capturing current view from screen...")
            
            # Directly save the axis scene to file (captures current view)
            # This avoids creating a new window/figure
            GLMakie.save(filename, ax_moment.blockscene, px_per_unit=2)
            
            println("✓ Moment space plot exported successfully!")
            println("  File size: $(round(filesize(filename)/1024, digits=1)) KB")
            println("="^70)
        catch e
            @warn "Moment space PNG export failed" exception=(e, catch_backtrace())
            println("Error: $e")
        end
    end
    
    # Isosurface controls with labels - wider sliders
    slider_iso1 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[1], width=200)
    slider_iso2 = slider_iso1  # Use same level
    slider_iso3 = slider_iso1  # Use same level
    slider_alpha = GLMakie.Slider(fig, range=0.3:0.1:1.0, startvalue=0.6, width=200)
    toggle_isosurface = GLMakie.Toggle(fig, active=true)
    toggle_streamlines = GLMakie.Toggle(fig, active=false)
    
    controls[5, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso Level", fontsize=9, halign=:left),
        slider_iso1;
        tellwidth=false, tellheight=false
    )
    controls[6, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Alpha", fontsize=9, halign=:left),
        slider_alpha;
        tellwidth=false, tellheight=false
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
                delete!(ax_physical, plot)
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
                p = GLMakie.contour!(ax_physical, x_lims, y_lims, z_lims, data,
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
                delete!(ax_physical, plot)
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
                p = GLMakie.streamplot!(ax_physical, 
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
    
    # Function to update moment space
    moment_plots = []
    slider_moment_threshold = nothing
    slider_boundary_alpha = nothing
    
    function update_moment_space!()
        if !has_std_moments || ax_moment === nothing
            return
        end
        
        # Clear previous plots
        for plot in moment_plots
            try
                delete!(ax_moment, plot)
            catch
            end
        end
        empty!(moment_plots)
        
        # Get current snapshot
        idx = time_slider.value[]
        S_field = snapshots[idx].S
        
        # Get correlation moments
        S110 = get_standardized_moment(S_field, "S110")
        S101 = get_standardized_moment(S_field, "S101")
        S011 = get_standardized_moment(S_field, "S011")
        
        # Compute correlation magnitude
        corr_mag = sqrt.(S110.^2 .+ S101.^2 .+ S011.^2)
        
        # Flatten arrays
        S110_flat = S110[:]
        S101_flat = S101[:]
        S011_flat = S011[:]
        corr_mag_flat = corr_mag[:]
        
        # Filter by threshold
        threshold = slider_moment_threshold !== nothing ? slider_moment_threshold.value[] : 0.01
        mask = corr_mag_flat .> threshold
        
        if sum(mask) > 0
            S110_filtered = S110_flat[mask]
            S101_filtered = S101_flat[mask]
            S011_filtered = S011_flat[mask]
            mag_filtered = corr_mag_flat[mask]
            
            # Draw point cloud in moment space
            p = GLMakie.scatter!(ax_moment, 
                               S110_filtered, S101_filtered, S011_filtered,
                               color=mag_filtered,
                               colormap=:viridis,
                               markersize=5,
                               alpha=0.6)
            push!(moment_plots, p)
        end
        
        # Draw coordinate axes
        p1 = GLMakie.lines!(ax_moment, [-1.0, 1.0], [0, 0], [0, 0], 
                     color=:red, linewidth=2, alpha=0.3)
        p2 = GLMakie.lines!(ax_moment, [0, 0], [-1.0, 1.0], [0, 0], 
                     color=:green, linewidth=2, alpha=0.3)
        p3 = GLMakie.lines!(ax_moment, [0, 0], [0, 0], [-1.0, 1.0], 
                     color=:blue, linewidth=2, alpha=0.3)
        push!(moment_plots, p1)
        push!(moment_plots, p2)
        push!(moment_plots, p3)
        
        # Draw |Δ₁| = 0 realizability boundary surface with proper transparency
        # Δ₁ = 1 + 2*S110*S101*S011 - S110² - S101² - S011² = 0
        # Use mesh/surface approach with observable for dynamic alpha control
        
        try
            # Get current boundary alpha from mapped slider (default 0.3)
            boundary_alpha_val = slider_boundary_alpha !== nothing ? slider_boundary_alpha[] : 0.3
            boundary_alpha_val = clamp(boundary_alpha_val, 0.0, 1.0)
            println("  Boundary alpha: $(round(boundary_alpha_val, digits=3)) ($(round(boundary_alpha_val*100, digits=1))%)")
            
            # Create 3D grid for isosurface
            n_grid = 60
            s_range = range(-1.0, 1.0, length=n_grid)
            
            # Compute Δ₁ at all grid points
            Delta1_volume = zeros(Float32, n_grid, n_grid, n_grid)
            @inbounds for (i, s110) in enumerate(s_range)
                @inbounds for (j, s101) in enumerate(s_range)
                    @inbounds for (k, s011) in enumerate(s_range)
                        Delta1_volume[i, j, k] = 1.0 + 2.0*s110*s101*s011 - s110^2 - s101^2 - s011^2
                    end
                end
            end
            
            # Create observable for boundary color with alpha
            boundary_color_obs = GLMakie.@lift(GLMakie.RGBAf(0.5, 0.5, 0.5, $slider_boundary_alpha))
            
            # Draw isosurface with observable color
            p_boundary = GLMakie.contour!(ax_moment,
                                         (-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0),
                                         Delta1_volume,
                                         levels=[0.0],
                                         color=boundary_color_obs,
                                         transparency=true,
                                         linewidth=0)
            
            push!(moment_plots, p_boundary)
            
            println("✓ Realizability boundary |Δ₁| = 0 displayed with dynamic transparency")
        catch e
            @warn "Could not compute realizability boundary" exception=e
            println("  Error: $e")
        end
        
        # Update moment space title
        ax_moment.title[] = latexstring("Moment Space - ", @sprintf("t=%.4f", snapshots[idx].t))
    end
    
    # Add moment threshold slider if we have standardized moments with label - wider
    if has_std_moments
        slider_moment_threshold = GLMakie.Slider(fig, range=0.001:0.001:0.5, startvalue=0.01, width=200)
        controls[7, 1] = GLMakie.vgrid!(
            GLMakie.Label(fig, "Min |S|", fontsize=9, halign=:left),
            slider_moment_threshold;
            tellwidth=false, tellheight=false
        )
        
        # Add realizability boundary opacity slider with fine control at low values
        # Use quadratic mapping: slider^2 gives more precision at low end
        slider_boundary_alpha_raw = GLMakie.Slider(fig, range=0.0:0.01:1.0, startvalue=0.55, width=200)
        controls[8, 1] = GLMakie.vgrid!(
            GLMakie.Label(fig, L"|Δ_1|~α", fontsize=9, halign=:left),
            slider_boundary_alpha_raw;
            tellwidth=false, tellheight=false
        )
        
        # Create mapped slider for fine low-end control
        slider_boundary_alpha = GLMakie.Observable(0.3)  # This will be the actual alpha value
        
        # Map raw slider value to alpha with quadratic scaling for fine low-end control
        GLMakie.on(slider_boundary_alpha_raw.value) do raw_val
            # Quadratic mapping: alpha = raw_val^2
            # This gives: 0.0->0.0, 0.1->0.01, 0.3->0.09, 0.5->0.25, 1.0->1.0
            mapped_alpha = raw_val^2
            slider_boundary_alpha[] = mapped_alpha
            update_moment_space!()
        end
        
        # Update moment space when threshold changes
        GLMakie.on(slider_moment_threshold.value) do val
            update_moment_space!()
        end
    end
    
    # Initial plots
    create_isosurfaces!()
    create_streamlines!()
    update_moment_space!()
    
    # Update plots when time slider changes
    GLMakie.on(time_slider.value) do val
        create_isosurfaces!()
        create_streamlines!()
        update_moment_space!()
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

