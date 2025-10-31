"""
Interactive 3D Time-Series Visualization with Streaming File Support

This viewer lazily loads snapshots from a JLD2 file on-demand as the user navigates,
rather than loading all snapshots into memory at once.

Note: This file is only loaded when visualization dependencies are available.
All required packages (GLMakie, FileIO, JLD2, LaTeXStrings, Dates, ColorSchemes)
are imported in the parent HyQMOM module.
"""

# Re-use imports from parent module (already imported conditionally)
using Printf
using ..LaTeXStrings  # Make @L_str macro available

# Import moment computation functions
import ..get_standardized_moment

"""
    interactive_3d_timeseries_streaming(filename, grid, params; kwargs...)

Launch an interactive 3D viewer that streams snapshots from a JLD2 file.

# Arguments
- `filename`: Path to JLD2 snapshot file
- `grid`: Grid structure with xm, ym, zm
- `params`: Simulation parameters

# Keyword Arguments
- `n_streamlines::Int=8`: Number of streamline seeds
- `vector_step::Int=4`: Subsampling for vector field
- `iso_levels::Vector{Float64}=[0.3, 0.5, 0.7]`: Isosurface levels
- `snapshot_mode::Symbol=:all`: Display mode (:all, :first, :last, :specific)
- `snapshot_number::Union{Int,Nothing}=nothing`: Specific snapshot number if mode is :specific
"""
function interactive_3d_timeseries_streaming(filename, grid, params;
                                             n_streamlines=8,
                                             vector_step=4,
                                             streamline_length=50,
                                             iso_levels=[0.3, 0.5, 0.7],
                                             snapshot_mode=:all,
                                             snapshot_number=nothing)
    
    println("\n" * "="^70)
    println("3D TIME-SERIES VIEWER (STREAMING MODE)")
    println("="^70)
    println("Features:")
    println("  * Snapshots loaded on-demand (low memory usage)")
    println("  * TRUE 3D isosurface contours")
    println("  * Velocity isosurfaces: Blue=positive, Red=negative")
    println("  * Time slider to navigate through evolution")
    println("="^70)
    
    # Open file and read metadata
    jld_file = JLD2.jldopen(filename, "r")
    n_snapshots = jld_file["meta/n_snapshots"]
    snap_keys = sort!(collect(keys(jld_file["snapshots"])))
    
    println("Loaded file: $filename")
    println("  $n_snapshots snapshots available")
    println("="^70)
    
    # Extract grid
    Nx = params.Nx
    Ny = params.Ny
    Nz = params.Nz
    xm = collect(grid.xm)
    ym = collect(grid.ym)
    zm = collect(grid.zm)
    
    # Check if snapshots have standardized moments
    first_snap = jld_file["snapshots/$(snap_keys[1])"]
    has_std_moments = haskey(first_snap, "S")
    
    println("\n" * "="^70)
    println("VIEWER LAYOUT")
    println("="^70)
    println("Has standardized moments (S field)? ", has_std_moments)
    
    # Determine initial snapshot index based on mode
    initial_snapshot_idx = if snapshot_mode == :first
        1
    elseif snapshot_mode == :last
        n_snapshots
    elseif snapshot_mode == :specific
        snapshot_number
    else
        1  # :all mode starts at first snapshot
    end
    
    # Create figure - always 3 columns for consistent layout
    println("Creating 3-column figure:")
    println("  Column 1: Physical space (x,y,z)")
    println("  Column 2: Moment space (S110, S101, S011)")
    println("  Column 3: Controls (fixed width)")
    # Compact window: 1500x600 total (1500 wide x 600 tall)
    fig = GLMakie.Figure(size=(1500, 600), fontsize=12,
                        fonts=(; regular="CMU Serif"))
    
    # Current snapshot index (observable)
    current_snapshot_idx = GLMakie.Observable(initial_snapshot_idx)
    
    # Left: Physical space (isosurfaces)
    ax_physical = GLMakie.Axis3(fig[1, 1], 
                                xlabel=L"x", ylabel=L"y", zlabel=L"z",
                                aspect=:data,
                                azimuth=0.3pi,
                                elevation=pi/8,
                                xticklabelsize=16, yticklabelsize=16, zticklabelsize=16,
                                xlabelsize=18, ylabelsize=18, zlabelsize=18)
    
    # Add coordinate axes to physical space (black, full opacity)
    x_center = (xm[1] + xm[end]) / 2
    y_center = (ym[1] + ym[end]) / 2
    z_center = (zm[1] + zm[end]) / 2
    GLMakie.lines!(ax_physical, [xm[1], xm[end]], [y_center, y_center], [z_center, z_center], 
                 color=:black, linewidth=2, alpha=1.0)
    GLMakie.lines!(ax_physical, [x_center, x_center], [ym[1], ym[end]], [z_center, z_center], 
                 color=:black, linewidth=2, alpha=1.0)
    GLMakie.lines!(ax_physical, [x_center, x_center], [y_center, y_center], [zm[1], zm[end]], 
                 color=:black, linewidth=2, alpha=1.0)
    
    # Middle: Moment space
    ax_moment = GLMakie.Axis3(fig[1, 2], 
                             xlabel=L"s_{110}", ylabel=L"s_{101}", zlabel=L"s_{011}",
                             aspect=:data,
                             azimuth=0.3pi,
                             elevation=pi/8,
                             limits=(-1, 1, -1, 1, -1, 1),
                             xticklabelsize=16, yticklabelsize=16, zticklabelsize=16,
                             xlabelsize=18, ylabelsize=18, zlabelsize=18)
    
    # Colorbar for moment space (below the 3D axis)
    # Fixed range from 0 to 1 (standardized moments are bounded)
    # Very light gray (almost white) → black for better contrast
    moment_cmap = [GLMakie.RGB(0.95, 0.95, 0.95), GLMakie.RGB(0.0, 0.0, 0.0)]
    cb_moment = GLMakie.Colorbar(fig[2, 2], 
                                 limits=(0.0, 1.0),
                                 colormap=moment_cmap,
                                 label=L"\Vert (s_{110},s_{101},s_{011}) \Vert_2",
                                 vertical=false,
                                 flipaxis=false,
                                 labelsize=18,
                                 ticklabelsize=16,
                                 height=18)
    
    println("="^70)
    
    # Control panel (far right)
    controls = fig[1, 3] = GLMakie.GridLayout(tellwidth=true)
    GLMakie.colsize!(fig.layout, 3, GLMakie.Fixed(250))
    
    # Set column gaps
    GLMakie.colgap!(fig.layout, 1, 10)
    GLMakie.colgap!(fig.layout, 2, 15)
    
    # Current quantity
    current_quantity = GLMakie.Observable("Density")
    
    # Quantity buttons
    btn_density = GLMakie.Button(fig, label="rho", fontsize=8)
    btn_u = GLMakie.Button(fig, label="u", fontsize=8)
    btn_v = GLMakie.Button(fig, label="v", fontsize=8)
    btn_w = GLMakie.Button(fig, label="w", fontsize=8)
    btn_pressure = GLMakie.Button(fig, label="P", fontsize=8)
    
    controls[1, 1] = GLMakie.hgrid!(btn_density, btn_u, btn_v, btn_w, btn_pressure; tellwidth=false)
    
    GLMakie.on(btn_density.clicks) do _
        current_quantity[] = "Density"
    end
    GLMakie.on(btn_u.clicks) do _
        current_quantity[] = "u velocity"
    end
    GLMakie.on(btn_v.clicks) do _
        current_quantity[] = "v velocity"
    end
    GLMakie.on(btn_w.clicks) do _
        current_quantity[] = "w velocity"
    end
    GLMakie.on(btn_pressure.clicks) do _
        current_quantity[] = "Pressure"
    end
    
    # Time slider and controls (only show if in :all mode)
    # Also create an observable for current snapshot index that works in both modes
    local time_slider
    if snapshot_mode == :all
        time_slider = GLMakie.Slider(fig, range=1:n_snapshots, startvalue=initial_snapshot_idx, width=200)
        btn_play = GLMakie.Button(fig, label=">", fontsize=8)
        btn_pause = GLMakie.Button(fig, label="||", fontsize=8)
        
        time_label = GLMakie.@lift(@sprintf("Snap %d/%d", $(time_slider.value), n_snapshots))
        
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
        
        # Use time_slider.value as the observable index
        current_snapshot_observable = time_slider.value
    else
        # Single snapshot mode - use a fixed observable at the initial index
        # Show which snapshot we're viewing
        if snapshot_mode == :first
            snap_label_text = "Showing: First snapshot"
        elseif snapshot_mode == :last
            snap_label_text = "Showing: Last snapshot"
        else
            snap_label_text = @sprintf("Showing: Snapshot %d/%d", initial_snapshot_idx, n_snapshots)
        end
        controls[2, 1] = GLMakie.Label(fig, snap_label_text, fontsize=10, halign=:center)
        
        is_playing = GLMakie.Observable(false)
        
        # Create an observable that stays fixed at the initial snapshot
        current_snapshot_observable = GLMakie.Observable(initial_snapshot_idx)
    end
    
    # Export button
    btn_export = GLMakie.Button(fig, label="Save PNG", fontsize=8)
    controls[7, 1] = btn_export
    
    GLMakie.on(btn_export.clicks) do _
        try
            idx = current_snapshot_observable[]
            snap_key = snap_keys[idx]
            snap = jld_file["snapshots/$snap_key"]
            t = snap["t"]
            timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
            quantity_short = replace(current_quantity[], " " => "_")
            
            filename_out = @sprintf("snapshot_%s_Nx%d_Ny%d_Nz%d_Kn%.2f_t%.4f_%s.png",
                               quantity_short, Nx, Ny, Nz, params.Kn, t, timestamp)
            
            println("\nExporting current view to: $filename_out")
            
            img = GLMakie.Makie.colorbuffer(fig.scene)
            FileIO.save(filename_out, img)
            
            println("[OK] Export complete!")
        catch e
            @error "Export failed" exception=(e, catch_backtrace())
        end
    end
    
    # Isosurface controls
    slider_iso1 = GLMakie.Slider(fig, range=0.1:0.05:0.9, startvalue=iso_levels[1], width=200)
    slider_alpha = GLMakie.Slider(fig, range=0.3:0.1:1.0, startvalue=0.6, width=200)
    
    # Function to compute quantities from M snapshot
    function compute_quantities(M_snapshot)
        rho = M_snapshot[:, :, :, 1]
        u = M_snapshot[:, :, :, 2] ./ rho
        v = M_snapshot[:, :, :, 6] ./ rho
        w = M_snapshot[:, :, :, 16] ./ rho
        
        C200 = M_snapshot[:, :, :, 3] ./ rho .- u.^2
        C020 = M_snapshot[:, :, :, 7] ./ rho .- v.^2
        C002 = M_snapshot[:, :, :, 17] ./ rho .- w.^2
        # Pressure: P = rho * (1/3 trace of velocity covariance)
        pressure = rho .* (C200 .+ C020 .+ C002) ./ 3.0
        
        return (rho=rho, u=u, v=v, w=w, pressure=pressure)
    end
    
    # Observable for current data (loads from file)
    current_data_obs = GLMakie.@lift begin
        idx = $(current_snapshot_observable)
        q = $(current_quantity)
        snap_key = snap_keys[idx]
        M = jld_file["snapshots/$snap_key/M"]
        quants = compute_quantities(M)
        
        if q == "Density"
            quants.rho
        elseif q == "u velocity"
            quants.u
        elseif q == "v velocity"
            quants.v
        elseif q == "w velocity"
            quants.w
        else # Pressure
            quants.pressure
        end
    end
    
    # Observable for isosurface value display (must come after current_data_obs)
    iso_value_text = GLMakie.@lift begin
        data = $(current_data_obs)
        q = $(current_quantity)
        iso_frac = $(slider_iso1.value)
        
        is_velocity = (q == "u velocity" || q == "v velocity" || q == "w velocity")
        data_min = minimum(data)
        data_max = maximum(data)
        data_absmax = maximum(abs.(data))
        
        # Get short label for quantity
        q_label = if q == "u velocity"
            "u"
        elseif q == "v velocity"
            "v"
        elseif q == "w velocity"
            "w"
        elseif q == "Density"
            "ρ"
        elseif q == "Pressure"
            "P"
        else
            q
        end
        
        if is_velocity
            # For velocities: show both positive and negative levels
            pos_level = iso_frac * data_absmax
            @sprintf("%s = ±%.4f", q_label, pos_level)
        else
            # For other quantities: show single level
            data_range = data_max - data_min
            level = data_min + iso_frac * data_range
            @sprintf("%s = %.4f", q_label, level)
        end
    end
    
    # Add isosurface controls to sidebar
    controls[4, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Iso Level", fontsize=9, halign=:left),
        slider_iso1,
        GLMakie.Label(fig, iso_value_text, fontsize=8, halign=:left, color=:gray);
        tellwidth=false
    )
    controls[5, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Alpha", fontsize=9, halign=:left),
        slider_alpha;
        tellwidth=false
    )
    
    # Storage for plots
    iso_plots = []
    
    # Prevent overlapping updates from sliders
    is_updating_isosurfaces = Ref(false)
    
    # Legend (persistent; we'll recreate it when needed but only once per update)
    current_legend = Ref{Union{Nothing, GLMakie.Legend}}(nothing)
    legend_layout = GLMakie.GridLayout(fig[2, 1])
    
    # Function to create isosurfaces
    function create_isosurfaces!()
        # Prevent concurrent updates
        if is_updating_isosurfaces[]
            return
        end
        is_updating_isosurfaces[] = true
        
        try
            for plot in iso_plots
                try
                    delete!(ax_physical, plot)
                catch
                end
            end
            empty!(iso_plots)
        
            data = current_data_obs[]
        q = current_quantity[]
        
        if any(isnan.(data)) || any(isinf.(data))
            @warn "Data contains NaN/Inf, skipping isosurfaces"
            return
        end
        
        is_velocity = (q == "u velocity" || q == "v velocity" || q == "w velocity")
        
        data_min = minimum(data)
        data_max = maximum(data)
        data_absmax = maximum(abs.(data))
        data_range = data_max - data_min
        
        if data_absmax < 1e-10 || data_range < 1e-10
            return
        end
        
        x_lims = (xm[1], xm[end])
        y_lims = (ym[1], ym[end])
        z_lims = (zm[1], zm[end])
        
        if is_velocity
            pos_level = slider_iso1.value[] * data_absmax
            neg_level = -slider_iso1.value[] * data_absmax
            
            levels = [pos_level, neg_level]
            # Dark blue for positive velocity, seagreen for negative velocity
            colors = [GLMakie.RGBf(0.0, 0.0, 0.545), :seagreen]
            alphas = [0.6, 0.6] .* slider_alpha.value[]
        else
            level = data_min + slider_iso1.value[] * data_range
            levels = [level]
            # Darkish blue for density and pressure
            colors = [GLMakie.RGBf(0.0, 0.0, 0.545)]
            alphas = [0.6] .* slider_alpha.value[]
        end
        
        # Build legend entries
        legend_entries = []
        
        for (idx, (level, color, alpha)) in enumerate(zip(levels, colors, alphas))
            if abs(level) < 1e-10
                continue
            end
            
            try
                p = GLMakie.contour!(ax_physical, x_lims, y_lims, z_lims, data,
                                    levels=[level],
                                    alpha=alpha,
                                    color=color)
                push!(iso_plots, p)
                
                # Format legend entry
                value_str = @sprintf("%.4f", level)
                if q == "Density"
                    entry = (color, L"\rho = %$(value_str)")
                elseif q == "u velocity"
                    entry = (color, L"u = %$(value_str)")
                elseif q == "v velocity"
                    entry = (color, L"v = %$(value_str)")
                elseif q == "w velocity"
                    entry = (color, L"w = %$(value_str)")
                elseif q == "Pressure"
                    entry = (color, L"P = %$(value_str)")
                else
                    entry = (color, @sprintf("Q = %.4f", level))
                end
                push!(legend_entries, entry)
            catch e
                if abs(level) > 1e-8
                    @warn "Contour failed at level $level" exception=(e,)
                end
            end
        end
        
            # Update legend - delete old one and create new one atomically
            if !isnothing(current_legend[])
                delete!(current_legend[])
            end
            
            if !isempty(legend_entries)
                legend_elements = [GLMakie.PolyElement(color=c, strokecolor=c, strokewidth=1) for (c, _) in legend_entries]
                legend_labels = [l for (_, l) in legend_entries]
                current_legend[] = GLMakie.Legend(legend_layout[1, 1], legend_elements, legend_labels,
                                                 orientation=:horizontal,
                                                 framevisible=false,
                                                 labelsize=14,
                                                 tellwidth=false,
                                                 tellheight=true,
                                                 patchsize=(30, 15))
            else
                current_legend[] = nothing
            end
        finally
            is_updating_isosurfaces[] = false
        end
    end
    
    # Function to update moment space
    moment_plots = []
    slider_moment_threshold = GLMakie.Slider(fig, range=0.001:0.001:0.5, startvalue=0.01, width=200)
    moment_range_text = GLMakie.Observable("|S| range: N/A")
    
    function update_moment_space!()
        if !has_std_moments
            return
        end
        
        for plot in moment_plots
            try
                delete!(ax_moment, plot)
            catch
            end
        end
        empty!(moment_plots)
        
        idx = current_snapshot_observable[]
        snap_key = snap_keys[idx]
        S_field = jld_file["snapshots/$snap_key/S"]
        
        S110 = get_standardized_moment(S_field, "S110")
        S101 = get_standardized_moment(S_field, "S101")
        S011 = get_standardized_moment(S_field, "S011")
        
        corr_mag = sqrt.(S110.^2 .+ S101.^2 .+ S011.^2)
        
        S110_flat = S110[:]
        S101_flat = S101[:]
        S011_flat = S011[:]
        corr_mag_flat = corr_mag[:]
        
        threshold = slider_moment_threshold.value[]
        mask = corr_mag_flat .> threshold
        
        if sum(mask) > 0
            S110_filtered = S110_flat[mask]
            S101_filtered = S101_flat[mask]
            S011_filtered = S011_flat[mask]
            mag_filtered = corr_mag_flat[mask]
            
            # Update range display
            min_mag = minimum(mag_filtered)
            max_mag = maximum(mag_filtered)
            moment_range_text[] = @sprintf("|S| range: %.3f - %.3f", min_mag, max_mag)
            
            # Low magnitude → very light gray, high magnitude → black
            # Use fixed colorrange from 0 to 1 for consistency with colorbar
            moment_scatter_cmap = [GLMakie.RGB(0.95, 0.95, 0.95), GLMakie.RGB(0.0, 0.0, 0.0)]
            p = GLMakie.scatter!(ax_moment, 
                               S110_filtered, S101_filtered, S011_filtered,
                               color=mag_filtered,
                               colormap=moment_scatter_cmap,
                               colorrange=(0.0, 1.0),
                               markersize=5,
                               alpha=0.6)
            push!(moment_plots, p)
        else
            moment_range_text[] = "|S| range: N/A"
        end
    end
    
    # Draw coordinate axes for moment space ONCE (outside update function)
    # Fixed black color with full opacity
    if has_std_moments
        GLMakie.lines!(ax_moment, [-1.0, 1.0], [0, 0], [0, 0], 
                     color=:black, linewidth=2, alpha=1.0)
        GLMakie.lines!(ax_moment, [0, 0], [-1.0, 1.0], [0, 0], 
                     color=:black, linewidth=2, alpha=1.0)
        GLMakie.lines!(ax_moment, [0, 0], [0, 0], [-1.0, 1.0], 
                     color=:black, linewidth=2, alpha=1.0)
    end
    
    # Draw |Delta_1| = 0 realizability boundary surface (transparent)
    # Delta_1 = 1 + 2*S110*S101*S011 - S110^2 - S101^2 - S011^2 = 0
    # This is the boundary of the realizable region in moment space
    # Draw this ONCE outside update function so slider doesn't affect it
    if has_std_moments
        try
            # Create a grid for the surface
            n_points = 50
            s1_range = range(-1, 1, length=n_points)
            s2_range = range(-1, 1, length=n_points)
            
            # We'll create the surface by solving for S011 given S110, S101
            # Rearranging: S011^2 - 2*S110*S101*S011 + (S110^2 + S101^2 - 1) = 0
            # Using quadratic formula: S011 = S110*S101 +/- sqrt((S110*S101)^2 - (S110^2 + S101^2 - 1))
            
            S110_grid = zeros(n_points, n_points)
            S101_grid = zeros(n_points, n_points)
            S011_grid_pos = zeros(n_points, n_points)
            S011_grid_neg = zeros(n_points, n_points)
            
            for (i, s110) in enumerate(s1_range)
                for (j, s101) in enumerate(s2_range)
                    S110_grid[i, j] = s110
                    S101_grid[i, j] = s101
                    
                    # Quadratic formula coefficients
                    # S011^2 - 2*a*b*S011 + (a^2 + b^2 - 1) = 0
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
            
            # Draw both sheets of the boundary surface with fixed low alpha
            GLMakie.surface!(ax_moment, 
                            S110_grid, S101_grid, S011_grid_pos,
                            color=:gray,
                            alpha=0.15,  # Fixed low transparency
                            transparency=true)
            
            GLMakie.surface!(ax_moment, 
                            S110_grid, S101_grid, S011_grid_neg,
                            color=:gray,
                            alpha=0.15,  # Fixed low transparency
                            transparency=true)
        catch e
            @warn "Could not compute realizability boundary" exception=e
        end
    end
    
    # Add moment threshold slider and range display
    controls[8, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Min |S|", fontsize=9, halign=:left),
        slider_moment_threshold,
        GLMakie.Label(fig, moment_range_text, fontsize=8, halign=:left, color=:gray);
        tellwidth=false
    )
    
    # Initial plots
    create_isosurfaces!()
    update_moment_space!()
    
    # Update plots when snapshot changes (only in :all mode with time slider)
    if snapshot_mode == :all
        GLMakie.on(current_snapshot_observable) do val
            create_isosurfaces!()
            update_moment_space!()
        end
    end
    
    # Update plots when quantity changes
    GLMakie.on(current_quantity) do q
        create_isosurfaces!()
    end
    
    # Update plots when sliders change (with error handling to prevent hanging)
    for slider in [slider_iso1, slider_alpha]
        GLMakie.on(slider.value) do val
            try
                create_isosurfaces!()
            catch e
                @warn "Isosurface update failed" exception=(e, catch_backtrace())
            end
        end
    end
    
    # Update moment space when threshold changes
    GLMakie.on(slider_moment_threshold.value) do val
        update_moment_space!()
    end
    
    # Animation loop for playback (only in :all mode)
    if snapshot_mode == :all
        GLMakie.on(is_playing) do playing
            if playing
                @async begin
                    while is_playing[] && current_snapshot_observable[] < n_snapshots
                        sleep(0.1)
                        if is_playing[]
                            current_snapshot_observable[] = current_snapshot_observable[] + 1
                        end
                    end
                    is_playing[] = false
                end
            end
        end
    end
    
    # Display the figure
    display(fig)
    
    println("\n" * "="^70)
    if snapshot_mode == :all
        println("TIME-SERIES VIEWER READY! (STREAMING MODE)")
        println("="^70)
        println("Layout:")
        println("  * Left: Physical space (x,y,z)")
        println("  * Middle: Moment space (S_1_1_0, S_1_0_1, S_0_1_1)")
        println("  * Right: Controls")
        println("\nControls:")
        println("  * Time slider steps through snapshots (loaded on-demand)")
        println("  * Click > Play to animate")
        println("  * Click quantity buttons to switch")
        println("  * Iso level/alpha sliders adjust appearance")
        println("  * Min |S| slider filters moment space")
    else
        println("SINGLE SNAPSHOT VIEWER READY!")
        println("="^70)
        if snapshot_mode == :first
            println("Viewing: First snapshot")
        elseif snapshot_mode == :last
            println("Viewing: Last snapshot")
        else
            println("Viewing: Snapshot $initial_snapshot_idx of $n_snapshots")
        end
        println("\nLayout:")
        println("  * Left: Physical space (x,y,z)")
        println("  * Middle: Moment space (S_1_1_0, S_1_0_1, S_0_1_1)")
        println("  * Right: Controls")
        println("\nControls:")
        println("  * Click quantity buttons to switch")
        println("  * Iso level/alpha sliders adjust appearance")
        println("  * Min |S| slider filters moment space")
    end
    println("\nPress Enter in terminal to close.")
    println("="^70)
    
    readline()
    
    # Close the JLD2 file
    close(jld_file)
    
    return fig
end

