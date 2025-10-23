"""
Interactive 3D Time-Series Visualization

This viewer allows stepping through simulation snapshots over time.
"""

import GLMakie
using Printf
using LaTeXStrings
using Dates

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
    if has_std_moments
        println("Creating 3-column figure (1800×1000):")
        println("  Column 1: Physical space (x,y,z)")
        println("  Column 2: Moment space (S110, S101, S011)")
        println("  Column 3: Controls")
        # Same width as standard, just split among 3 columns, minimal gaps
        # Use LaTeX fonts for all plot text
        fig = GLMakie.Figure(size=(1600, 700), fontsize=12,
                            fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
        # Reduce column gaps to maximize plot space
        GLMakie.colgap!(fig.layout, 5)  # 5 pixels between columns
    else
        println("Creating 2-column figure (1800×1000)")
        # Default size: 1800x1000
        # Use LaTeX fonts for all plot text
        fig = GLMakie.Figure(size=(1600, 700), fontsize=12,
                            fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
        GLMakie.colgap!(fig.layout, 5)  # 5 pixels between columns
    end
    
    # Current snapshot index (observable) - needed for moment space title
    current_snapshot_idx = GLMakie.Observable(1)
    
    # Left: Physical space (isosurfaces)
    ax_physical = GLMakie.Axis3(fig[1, 1], 
                                xlabel=L"x", ylabel=L"y", zlabel=L"z",
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
    btn_density = GLMakie.Button(fig, label="ρ", fontsize=8)
    btn_u = GLMakie.Button(fig, label="U", fontsize=8)
    btn_v = GLMakie.Button(fig, label="V", fontsize=8)
    btn_w = GLMakie.Button(fig, label="W", fontsize=8)
    btn_pressure = GLMakie.Button(fig, label="P", fontsize=8)
    
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
    btn_play = GLMakie.Button(fig, label=">", fontsize=8)
    btn_pause = GLMakie.Button(fig, label="||", fontsize=8)
    
    time_label = GLMakie.@lift(@sprintf("Snap %d/%d", $(time_slider.value), length(snapshots)))
    
    controls[2, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, time_label, fontsize=9, halign=:left),
        time_slider;
        tellwidth=false
    )
    controls[3, 1] = GLMakie.hgrid!(btn_play, btn_pause; tellwidth=false)
    
    is_playing = GLMakie.Observable(false)
    
    GLMakie.on(btn_play.clicks) do _
        is_playing[] = true
    end
    GLMakie.on(btn_pause.clicks) do _
        is_playing[] = false
    end
    
    # Export button - save current figure at high resolution
    btn_export = GLMakie.Button(fig, label="Save PNG", fontsize=8)
    controls[7, 1] = btn_export
    
    GLMakie.on(btn_export.clicks) do _
        try
            idx = time_slider.value[]
            snap = snapshots[idx]
            timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
            quantity_short = replace(current_quantity[], " " => "_")
            filename = @sprintf("snapshot_%s_t%.4f_%s.png", quantity_short, snap.t, timestamp)
            
            println("\nExporting current view to: $filename")
            println("  Resolution: $(fig.scene.viewport[].widths .* 8) pixels")
            
            # Save at 4× native resolution for publication quality
            # This gives approximately 300 DPI for a figure that's ~5 inches wide
            GLMakie.save(filename, fig; px_per_unit=8)
            
            println("✓ Export complete!")
            println("  File saved: $filename")
        catch e
            @error "Export failed" exception=(e, catch_backtrace())
            println("Failed to export figure.")
        end
    end
    
    # Isosurface controls with labels - wider sliders
    slider_iso1 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[1], width=200)
    slider_iso2 = slider_iso1  # Use same level
    slider_iso3 = slider_iso1  # Use same level
    slider_alpha = GLMakie.Slider(fig, range=0.3:0.1:1.0, startvalue=0.6, width=200)
    toggle_isosurface = GLMakie.Toggle(fig, active=true)
    toggle_streamlines = GLMakie.Toggle(fig, active=false)
    
    controls[4, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso Level", fontsize=9, halign=:left),
        slider_iso1;
        tellwidth=false
    )
    controls[5, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Alpha", fontsize=9, halign=:left),
        slider_alpha;
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
        
        # Draw |Δ₁| = 0 realizability boundary surface (transparent)
        # Δ₁ = 1 + 2*S110*S101*S011 - S110² - S101² - S011² = 0
        # This is the boundary of the realizable region in moment space
        
        try
            # Create a grid for the surface
            n_points = 50
            s1_range = range(-1, 1, length=n_points)
            s2_range = range(-1, 1, length=n_points)
            
            # We'll create the surface by solving for S011 given S110, S101
            # Rearranging: S011² - 2*S110*S101*S011 + (S110² + S101² - 1) = 0
            # Using quadratic formula: S011 = S110*S101 ± sqrt((S110*S101)² - (S110² + S101² - 1))
            
            S110_grid = zeros(n_points, n_points)
            S101_grid = zeros(n_points, n_points)
            S011_grid_pos = zeros(n_points, n_points)
            S011_grid_neg = zeros(n_points, n_points)
            
            for (i, s110) in enumerate(s1_range)
                for (j, s101) in enumerate(s2_range)
                    S110_grid[i, j] = s110
                    S101_grid[i, j] = s101
                    
                    # Quadratic formula coefficients
                    # S011² - 2*a*b*S011 + (a² + b² - 1) = 0
                    discriminant = (s110 * s101)^2 - (s110^2 + s101^2 - 1)
                    
                    if discriminant >= 0
                        sqrt_disc = sqrt(discriminant)
                        S011_grid_pos[i, j] = s110 * s101 + sqrt_disc
                        S011_grid_neg[i, j] = s110 * s101 - sqrt_disc
                    else
                        # No real solution - mark as NaN (won't plot)
                        S011_grid_pos[i, j] = NaN
                        S011_grid_neg[i, j] = NaN
                    end
                end
            end
            
            # Clamp to [-1, 1] range
            S011_grid_pos = clamp.(S011_grid_pos, -1, 1)
            S011_grid_neg = clamp.(S011_grid_neg, -1, 1)
            
            # Draw both sheets of the boundary surface
            p_boundary_pos = GLMakie.surface!(ax_moment, 
                                            S110_grid, S101_grid, S011_grid_pos,
                                            color=:gray,
                                            alpha=0.15,  # Very transparent
                                            transparency=true)
            
            p_boundary_neg = GLMakie.surface!(ax_moment, 
                                            S110_grid, S101_grid, S011_grid_neg,
                                            color=:gray,
                                            alpha=0.15,  # Very transparent
                                            transparency=true)
            
            push!(moment_plots, p_boundary_pos)
            push!(moment_plots, p_boundary_neg)
            
            println("✓ Realizability boundary |Δ₁| = 0 displayed")
        catch e
            @warn "Could not compute realizability boundary" exception=e
        end
        
        # Title removed for cleaner visualization
    end
    
    # Add moment threshold slider if we have standardized moments with label - wider
    if has_std_moments
        slider_moment_threshold = GLMakie.Slider(fig, range=0.001:0.001:0.5, startvalue=0.01, width=200)
        controls[6, 1] = GLMakie.vgrid!(
            GLMakie.Label(fig, "Min |S|", fontsize=9, halign=:left),
            slider_moment_threshold;
            tellwidth=false
        )
        
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
    println("  • Mouse: drag to rotate, scroll to zoom camera")
    println("  • Resize window to zoom entire plot (including axes/labels)")
    println("  • Click 'Save PNG' button to export current view (high resolution)")
    println("\nPress Enter in terminal to close.")
    println("="^70)
    
    readline()
    
    return fig
end

