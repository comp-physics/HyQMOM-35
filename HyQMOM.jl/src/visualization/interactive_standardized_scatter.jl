"""
Interactive 3D Scatterplot Visualization for Standardized Moments

This viewer displays standardized moments as 3D scatterplots/point clouds,
useful for visualizing correlations and moment distributions.
"""

import GLMakie
using Printf
using LaTeXStrings

"""
    interactive_standardized_scatter(snapshot, grid; kwargs...)

Launch an interactive 3D scatterplot viewer for standardized moments.

# Arguments
- `snapshot`: NamedTuple with fields `M`, `S` (standardized moments), `t`, `step`
- `grid`: Grid structure with `xm`, `ym`, `zm`

# Keyword Arguments
- `threshold::Float64=0.1`: Only show points where |moment| > threshold
- `subsample::Int=1`: Subsample grid points (1=all, 2=every other, etc.)
- `markersize::Float64=3.0`: Size of scatter points
- `colormap::Symbol=:RdBu`: Colormap for scatter points

# Example
```julia
using HyQMOM, JLD2

# Load snapshots with standardized moments
@load "snapshots.jld2" snapshots grid

# View standardized moments at snapshot 5
interactive_standardized_scatter(snapshots[5], grid)
```

The viewer allows you to:
- Select which standardized moment to display (S110, S101, S022, etc.)
- View as 3D scatterplot colored by moment value
- Filter by threshold to see only significant values
- Subsample for better performance with large grids
"""
function interactive_standardized_scatter(snapshot, grid;
                                         threshold=0.1,
                                         subsample=1,
                                         markersize=3.0,
                                         colormap=:RdBu)
    
    # Check if standardized moments are available
    if !haskey(snapshot, :S)
        error("Snapshot does not contain standardized moments (S field). " *
              "Run simulation with save_standardized_moments=true")
    end
    
    println("\n" * "="^70)
    println("3D SCATTERPLOT VIEWER - STANDARDIZED MOMENTS")
    println("="^70)
    println("Features:")
    println("  • View standardized moments as 3D point clouds")
    println("  • Color by moment value")
    println("  • Filter by threshold")
    println("  • Explore moment correlations and distributions")
    println("="^70)
    
    S_field = snapshot.S
    Nx, Ny, Nz, Nmom = size(S_field)
    
    xm = collect(grid.xm)
    ym = collect(grid.ym)
    zm = collect(grid.zm)
    
    # Get full domain extents
    x_min, x_max = minimum(xm), maximum(xm)
    y_min, y_max = minimum(ym), maximum(ym)
    z_min, z_max = minimum(zm), maximum(zm)
    
    # Create figure with LaTeX font rendering
    fig = GLMakie.Figure(size=(1800, 1000), fontsize=12,
                        fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
    
    # Main 3D axis with fixed limits for full domain
    ax = GLMakie.Axis3(fig[1:3, 1:2], 
                       xlabel=L"x", ylabel=L"y", zlabel=L"z",
                       aspect=:data,
                       azimuth=0.3π,
                       elevation=π/8,
                       limits=(x_min, x_max, y_min, y_max, z_min, z_max),
                       xticklabelsize=11, yticklabelsize=11, zticklabelsize=11,
                       xlabelsize=13, ylabelsize=13, zlabelsize=13)
    
    # Control panel
    controls = fig[1:3, 3] = GLMakie.GridLayout()
    
    # Current moment selection
    current_moment = GLMakie.Observable("S110")
    
    # Moment selector buttons
    controls[1, 1] = GLMakie.Label(fig, "Standardized Moments:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    # Group moments by order
    controls[2, 1] = GLMakie.Label(fig, "2nd Order (Correlations):", fontsize=12, 
                                   halign=:left, tellwidth=false)
    
    btn_s110 = GLMakie.Button(fig, label=L"S_{110}", fontsize=10, width=90)
    btn_s101 = GLMakie.Button(fig, label=L"S_{101}", fontsize=10, width=90)
    btn_s011 = GLMakie.Button(fig, label=L"S_{011}", fontsize=10, width=90)
    controls[3, 1] = GLMakie.hgrid!(btn_s110, btn_s101, btn_s011; tellwidth=false)
    
    controls[4, 1] = GLMakie.Label(fig, "4th Order (Anisotropy):", fontsize=12, 
                                   halign=:left, tellwidth=false)
    
    btn_s220 = GLMakie.Button(fig, label=L"S_{220}", fontsize=10, width=90)
    btn_s202 = GLMakie.Button(fig, label=L"S_{202}", fontsize=10, width=90)
    btn_s022 = GLMakie.Button(fig, label=L"S_{022}", fontsize=10, width=90)
    controls[5, 1] = GLMakie.hgrid!(btn_s220, btn_s202, btn_s022; tellwidth=false)
    
    controls[6, 1] = GLMakie.Label(fig, "3rd Order (Skewness):", fontsize=12, 
                                   halign=:left, tellwidth=false)
    
    btn_s300 = GLMakie.Button(fig, label=L"S_{300}", fontsize=10, width=90)
    btn_s030 = GLMakie.Button(fig, label=L"S_{030}", fontsize=10, width=90)
    btn_s003 = GLMakie.Button(fig, label=L"S_{003}", fontsize=10, width=90)
    controls[7, 1] = GLMakie.hgrid!(btn_s300, btn_s030, btn_s003; tellwidth=false)
    
    btn_s210 = GLMakie.Button(fig, label=L"S_{210}", fontsize=10, width=60)
    btn_s201 = GLMakie.Button(fig, label=L"S_{201}", fontsize=10, width=60)
    btn_s120 = GLMakie.Button(fig, label=L"S_{120}", fontsize=10, width=60)
    btn_s102 = GLMakie.Button(fig, label=L"S_{102}", fontsize=10, width=60)
    controls[8, 1] = GLMakie.hgrid!(btn_s210, btn_s201, btn_s120, btn_s102; tellwidth=false)
    
    btn_s021 = GLMakie.Button(fig, label=L"S_{021}", fontsize=10, width=60)
    btn_s012 = GLMakie.Button(fig, label=L"S_{012}", fontsize=10, width=60)
    btn_s111 = GLMakie.Button(fig, label=L"S_{111}", fontsize=10, width=60)
    controls[9, 1] = GLMakie.hgrid!(btn_s021, btn_s012, btn_s111; tellwidth=false)
    
    controls[10, 1] = GLMakie.Label(fig, "4th Order (Kurtosis):", fontsize=12, 
                                    halign=:left, tellwidth=false)
    
    btn_s400 = GLMakie.Button(fig, label=L"S_{400}", fontsize=10, width=90)
    btn_s040 = GLMakie.Button(fig, label=L"S_{040}", fontsize=10, width=90)
    btn_s004 = GLMakie.Button(fig, label=L"S_{004}", fontsize=10, width=90)
    controls[11, 1] = GLMakie.hgrid!(btn_s400, btn_s040, btn_s004; tellwidth=false)
    
    # Connect buttons to observable
    for (btn, name) in [(btn_s110, "S110"), (btn_s101, "S101"), (btn_s011, "S011"),
                        (btn_s220, "S220"), (btn_s202, "S202"), (btn_s022, "S022"),
                        (btn_s300, "S300"), (btn_s030, "S030"), (btn_s003, "S003"),
                        (btn_s210, "S210"), (btn_s201, "S201"), (btn_s120, "S120"),
                        (btn_s102, "S102"), (btn_s021, "S021"), (btn_s012, "S012"),
                        (btn_s111, "S111"), (btn_s400, "S400"), (btn_s040, "S040"),
                        (btn_s004, "S004")]
        GLMakie.on(btn.clicks) do n
            current_moment[] = name
            println("Switched to: $name")
        end
    end
    
    # Visualization controls
    controls[13, 1] = GLMakie.Label(fig, "Display Options:", fontsize=14, 
                                    halign=:left, font=:bold, tellwidth=false)
    
    # Threshold slider
    slider_threshold = GLMakie.Slider(fig, range=0.0:0.05:1.0, startvalue=threshold)
    controls[14, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, @sprintf("Threshold: %.2f", threshold), fontsize=12, halign=:left),
        slider_threshold;
        tellwidth=false
    )
    
    # Marker size slider
    slider_markersize = GLMakie.Slider(fig, range=1.0:0.5:10.0, startvalue=markersize)
    controls[15, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, @sprintf("Point size: %.1f", markersize), fontsize=12, halign=:left),
        slider_markersize;
        tellwidth=false
    )
    
    # Subsample slider
    slider_subsample = GLMakie.Slider(fig, range=1:1:5, startvalue=subsample)
    controls[16, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, @sprintf("Subsample: %d", subsample), fontsize=12, halign=:left),
        slider_subsample;
        tellwidth=false
    )
    
    # Toggle for showing only positive/negative/both
    toggle_positive = GLMakie.Toggle(fig, active=true)
    toggle_negative = GLMakie.Toggle(fig, active=true)
    controls[17, 1] = GLMakie.hgrid!(
        toggle_positive, GLMakie.Label(fig, "Positive", fontsize=11, halign=:left),
        toggle_negative, GLMakie.Label(fig, "Negative", fontsize=11, halign=:left);
        tellwidth=false
    )
    
    # Update function with proper cleanup and update guard
    updating = Ref(false)  # Prevent recursive updates
    
    function update_scatter!()
        # Guard against recursive calls
        if updating[]
            return
        end
        
        try
            updating[] = true
            
            # Clear all plots from axis to prevent buildup
            GLMakie.empty!(ax)
            
            # Get current moment field
            moment_name = current_moment[]
            moment_data = get_standardized_moment(S_field, moment_name)
        
        # Subsample
        step = Int(slider_subsample.value[])
        indices_x = 1:step:Nx
        indices_y = 1:step:Ny
        indices_z = 1:step:Nz
        
        moment_sub = moment_data[indices_x, indices_y, indices_z]
        xm_sub = xm[indices_x]
        ym_sub = ym[indices_y]
        zm_sub = zm[indices_z]
        
        # Create 3D grid of points
        xs = [x for x in xm_sub, y in ym_sub, z in zm_sub][:]
        ys = [y for x in xm_sub, y in ym_sub, z in zm_sub][:]
        zs = [z for x in xm_sub, y in ym_sub, z in zm_sub][:]
        values = moment_sub[:]
        
        # Apply threshold and sign filters
        thresh = slider_threshold.value[]
        show_pos = toggle_positive.active[]
        show_neg = toggle_negative.active[]
        
        mask = abs.(values) .> thresh
        if !show_pos
            mask = mask .&& (values .<= 0)
        end
        if !show_neg
            mask = mask .&& (values .>= 0)
        end
        
        n_points = sum(mask)
        
        if n_points == 0
            println("No points meet threshold criteria")
            return
        end
        
        xs_filtered = xs[mask]
        ys_filtered = ys[mask]
        zs_filtered = zs[mask]
        values_filtered = values[mask]
        
        # Create scatter plot
        marker_sz = slider_markersize.value[]
        
            GLMakie.scatter!(ax,
                            xs_filtered, ys_filtered, zs_filtered,
                            color=values_filtered,
                            colormap=colormap,
                            markersize=marker_sz,
                            alpha=0.6)
            
            println("Displayed $n_points points ($(round(100*n_points/length(values), digits=1))% of total)")
        finally
            updating[] = false
        end
    end
    
    # Initial scatter
    update_scatter!()
    
    # Connect updates
    GLMakie.on(current_moment) do _
        update_scatter!()
    end
    
    GLMakie.on(slider_threshold.value) do _
        update_scatter!()
    end
    
    GLMakie.on(slider_markersize.value) do _
        update_scatter!()
    end
    
    GLMakie.on(slider_subsample.value) do _
        update_scatter!()
    end
    
    GLMakie.on(toggle_positive.active) do _
        update_scatter!()
    end
    
    GLMakie.on(toggle_negative.active) do _
        update_scatter!()
    end
    
    # Colorbar
    cb = GLMakie.Colorbar(fig[1:3, 4],
                         label="Moment Value",
                         colormap=colormap,
                         limits=(-1.0, 1.0),
                         width=25)
    
    # Statistics panel
    stats_text = GLMakie.@lift(begin
        moment_data = get_standardized_moment(S_field, $current_moment)
        name = $current_moment
        thresh = slider_threshold.value[]
        sub = slider_subsample.value[]
        
        n_above_thresh = sum(abs.(moment_data) .> thresh)
        pct_above = 100 * n_above_thresh / length(moment_data)
        
        @sprintf("""
        %s Statistics:
          Min:    %+.4f
          Max:    %+.4f
          Mean:   %+.4f
          Std:    %.4f
          
        Filtering:
          |%s| > %.2f: %d points (%.1f%%)
          Subsample: every %d points
          
        Grid: %dx%dx%d
        Time: t = %.4f (step %d)
        """, name, minimum(moment_data), maximum(moment_data),
             sum(moment_data)/length(moment_data), 
             sqrt(sum((moment_data .- sum(moment_data)/length(moment_data)).^2) / length(moment_data)),
             name, thresh, n_above_thresh, pct_above, sub,
             Nx, Ny, Nz, snapshot.t, snapshot.step)
    end)
    
    GLMakie.Label(fig[4, 1:4], stats_text, fontsize=11, halign=:left,
                 padding=(10, 10, 10, 10), tellwidth=false)
    
    println("\n" * "="^70)
    println("SCATTERPLOT VIEWER READY!")
    println("="^70)
    println("\nControls:")
    println("  • Click moment buttons to switch quantities")
    println("  • Adjust threshold to filter small values")
    println("  • Adjust point size for visibility")
    println("  • Subsample for better performance")
    println("  • Toggle positive/negative to filter by sign")
    println("  • Mouse: drag to rotate, scroll to zoom")
    println("\nPhysical Interpretation:")
    println("  S110: u-v velocity correlation (xy shear)")
    println("  S101: u-w velocity correlation (xz shear)")
    println("  S011: v-w velocity correlation (yz shear)")
    println("  S220, S202, S022: Temperature/energy anisotropy")
    println("  S300, S030, S003: Directional skewness")
    println("  S400, S040, S004: Directional kurtosis (Gaussian=3)")
    println("\nPress Enter to close.")
    println("="^70)
    
    display(fig)
    readline()
    
    return fig
end

"""
    interactive_standardized_scatter(snapshots::Vector, grid; kwargs...)

Launch an interactive 3D scatterplot viewer with time slider for standardized moments.

# Arguments
- `snapshots`: Vector of NamedTuples with fields `M`, `S` (standardized moments), `t`, `step`
- `grid`: Grid structure with `xm`, `ym`, `zm`

# Keyword Arguments
- `threshold::Float64=0.1`: Only show points where |moment| > threshold
- `subsample::Int=1`: Subsample grid points (1=all, 2=every other, etc.)
- `markersize::Float64=3.0`: Size of scatter points
- `colormap::Symbol=:RdBu`: Colormap for scatter points

# Example
```julia
using HyQMOM, JLD2

# Load snapshots with standardized moments
@load "snapshots.jld2" snapshots grid

# View all snapshots with time slider
interactive_standardized_scatter(snapshots, grid)
```

The viewer allows you to:
- Step through time using the time slider
- Select which standardized moment to display (S110, S101, S022, etc.)
- View as 3D scatterplot colored by moment value
- Filter by threshold to see only significant values
- Subsample for better performance with large grids
"""
function interactive_standardized_scatter(snapshots::Vector, grid;
                                         threshold=0.1,
                                         subsample=1,
                                         markersize=3.0,
                                         colormap=:RdBu)
    
    # Check if standardized moments are available
    if !haskey(snapshots[1], :S)
        error("Snapshots do not contain standardized moments (S field). " *
              "Run simulation with save_standardized_moments=true")
    end
    
    println("\n" * "="^70)
    println("3D SCATTERPLOT VIEWER - STANDARDIZED MOMENTS (TIME SERIES)")
    println("="^70)
    println("Features:")
    println("  • Step through simulation snapshots over time")
    println("  • View standardized moments as 3D point clouds")
    println("  • Color by moment value")
    println("  • Filter by threshold")
    println("  • Explore moment correlations and distributions")
    println("="^70)
    println("Loaded $(length(snapshots)) snapshots")
    println("  Time range: $(snapshots[1].t) to $(snapshots[end].t)")
    println("="^70)
    
    # Get grid info
    xm = collect(grid.xm)
    ym = collect(grid.ym)
    zm = collect(grid.zm)
    Nx, Ny, Nz = length(xm), length(ym), length(zm)
    
    # Get full domain extents
    x_min, x_max = minimum(xm), maximum(xm)
    y_min, y_max = minimum(ym), maximum(ym)
    z_min, z_max = minimum(zm), maximum(zm)
    
    # Create figure with LaTeX font rendering
    fig = GLMakie.Figure(size=(1800, 1000), fontsize=12,
                        fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
    
    # Current snapshot index
    current_snapshot_idx = GLMakie.Observable(1)
    
    # Main 3D axis with fixed limits for full domain
    ax = GLMakie.Axis3(fig[1:3, 1:2], 
                       xlabel=L"x", ylabel=L"y", zlabel=L"z",
                       aspect=:data,
                       azimuth=0.3π,
                       elevation=π/8,
                       limits=(x_min, x_max, y_min, y_max, z_min, z_max),
                       xticklabelsize=11, yticklabelsize=11, zticklabelsize=11,
                       xlabelsize=13, ylabelsize=13, zlabelsize=13)
    
    # Control panel
    controls = fig[1:3, 3] = GLMakie.GridLayout()
    
    # Time slider at the top
    controls[1, 1] = GLMakie.Label(fig, "Time Control:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    time_slider = GLMakie.Slider(fig, range=1:length(snapshots), startvalue=1)
    time_label = GLMakie.Label(fig, GLMakie.@lift(@sprintf("Snapshot %d/%d, t=%.4f", 
                                                          $(time_slider.value)[], 
                                                          length(snapshots),
                                                          snapshots[$(time_slider.value)[]].t)),
                              fontsize=11, halign=:left)
    controls[2, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Time:", fontsize=12, halign=:left),
        time_slider,
        time_label;
        tellwidth=false
    )
    
    # Update current snapshot when slider changes
    GLMakie.on(time_slider.value) do idx
        current_snapshot_idx[] = Int(idx)
    end
    
    # Current moment selection
    current_moment = GLMakie.Observable("S110")
    
    # Moment selector buttons
    controls[3, 1] = GLMakie.Label(fig, "Standardized Moments:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    # Group moments by order
    controls[4, 1] = GLMakie.Label(fig, "2nd Order (Correlations):", fontsize=12, 
                                   halign=:left, tellwidth=false)
    
    btn_s110 = GLMakie.Button(fig, label=L"S_{110}", fontsize=10, width=90)
    btn_s101 = GLMakie.Button(fig, label=L"S_{101}", fontsize=10, width=90)
    btn_s011 = GLMakie.Button(fig, label=L"S_{011}", fontsize=10, width=90)
    controls[5, 1] = GLMakie.hgrid!(btn_s110, btn_s101, btn_s011; tellwidth=false)
    
    controls[6, 1] = GLMakie.Label(fig, "4th Order (Anisotropy):", fontsize=12, 
                                   halign=:left, tellwidth=false)
    
    btn_s220 = GLMakie.Button(fig, label=L"S_{220}", fontsize=10, width=90)
    btn_s202 = GLMakie.Button(fig, label=L"S_{202}", fontsize=10, width=90)
    btn_s022 = GLMakie.Button(fig, label=L"S_{022}", fontsize=10, width=90)
    controls[7, 1] = GLMakie.hgrid!(btn_s220, btn_s202, btn_s022; tellwidth=false)
    
    controls[8, 1] = GLMakie.Label(fig, "3rd Order (Skewness):", fontsize=12, 
                                   halign=:left, tellwidth=false)
    
    btn_s300 = GLMakie.Button(fig, label=L"S_{300}", fontsize=10, width=90)
    btn_s030 = GLMakie.Button(fig, label=L"S_{030}", fontsize=10, width=90)
    btn_s003 = GLMakie.Button(fig, label=L"S_{003}", fontsize=10, width=90)
    controls[9, 1] = GLMakie.hgrid!(btn_s300, btn_s030, btn_s003; tellwidth=false)
    
    controls[10, 1] = GLMakie.Label(fig, "4th Order (Kurtosis):", fontsize=12, 
                                    halign=:left, tellwidth=false)
    
    btn_s400 = GLMakie.Button(fig, label=L"S_{400}", fontsize=10, width=90)
    btn_s040 = GLMakie.Button(fig, label=L"S_{040}", fontsize=10, width=90)
    btn_s004 = GLMakie.Button(fig, label=L"S_{004}", fontsize=10, width=90)
    controls[11, 1] = GLMakie.hgrid!(btn_s400, btn_s040, btn_s004; tellwidth=false)
    
    # Connect buttons to moment selection
    moment_buttons = [
        ("S110", btn_s110), ("S101", btn_s101), ("S011", btn_s011),
        ("S220", btn_s220), ("S202", btn_s202), ("S022", btn_s022),
        ("S300", btn_s300), ("S030", btn_s030), ("S003", btn_s003),
        ("S400", btn_s400), ("S040", btn_s040), ("S004", btn_s004)
    ]
    
    for (name, btn) in moment_buttons
        GLMakie.on(btn.clicks) do n
            current_moment[] = name
            println("Switched to: $name")
        end
    end
    
    # Visualization controls
    controls[12, 1] = GLMakie.Label(fig, "Display Options:", fontsize=14, 
                                    halign=:left, font=:bold, tellwidth=false)
    
    # Threshold slider
    slider_threshold = GLMakie.Slider(fig, range=0.0:0.05:1.0, startvalue=threshold)
    controls[13, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, @sprintf("Threshold: %.2f", threshold), fontsize=12, halign=:left),
        slider_threshold;
        tellwidth=false
    )
    
    # Marker size slider
    slider_markersize = GLMakie.Slider(fig, range=1.0:0.5:10.0, startvalue=markersize)
    controls[14, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, @sprintf("Point size: %.1f", markersize), fontsize=12, halign=:left),
        slider_markersize;
        tellwidth=false
    )
    
    # Subsample slider
    slider_subsample = GLMakie.Slider(fig, range=1:1:5, startvalue=subsample)
    controls[15, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, @sprintf("Subsample: %d", subsample), fontsize=12, halign=:left),
        slider_subsample;
        tellwidth=false
    )
    
    # Toggle for showing only positive/negative/both
    toggle_positive = GLMakie.Toggle(fig, active=true)
    toggle_negative = GLMakie.Toggle(fig, active=true)
    controls[16, 1] = GLMakie.hgrid!(
        GLMakie.Label(fig, "Show:", fontsize=12),
        GLMakie.Label(fig, "Positive", fontsize=11),
        toggle_positive,
        GLMakie.Label(fig, "Negative", fontsize=11),
        toggle_negative;
        tellwidth=false
    )
    
    # Update function with proper cleanup and update guard
    updating = Ref(false)  # Prevent recursive updates
    
    function update_scatter!()
        # Guard against recursive calls
        if updating[]
            return
        end
        
        try
            updating[] = true
            
            # Get current snapshot
            snap_idx = current_snapshot_idx[]
            S_field = snapshots[snap_idx].S
            
            # Clear all plots from axis to prevent buildup
            GLMakie.empty!(ax)
            
            # Get current moment field
            moment_name = current_moment[]
            moment_data = get_standardized_moment(S_field, moment_name)
            
            # Subsample
            step = Int(slider_subsample.value[])
            indices_x = 1:step:Nx
            indices_y = 1:step:Ny
            indices_z = 1:step:Nz
            
            moment_sub = moment_data[indices_x, indices_y, indices_z]
            xm_sub = xm[indices_x]
            ym_sub = ym[indices_y]
            zm_sub = zm[indices_z]
            
            # Create 3D grid of points
            xs = [x for x in xm_sub, y in ym_sub, z in zm_sub][:]
            ys = [y for x in xm_sub, y in ym_sub, z in zm_sub][:]
            zs = [z for x in xm_sub, y in ym_sub, z in zm_sub][:]
            values = moment_sub[:]
            
            # Apply threshold and sign filters
            thresh = slider_threshold.value[]
            show_pos = toggle_positive.active[]
            show_neg = toggle_negative.active[]
            
            mask = abs.(values) .> thresh
            if !show_pos
                mask = mask .&& (values .<= 0)
            end
            if !show_neg
                mask = mask .&& (values .>= 0)
            end
            
            n_points = sum(mask)
            
            if n_points == 0
                println("No points meet threshold criteria")
                return
            end
            
            xs_filtered = xs[mask]
            ys_filtered = ys[mask]
            zs_filtered = zs[mask]
            values_filtered = values[mask]
            
            # Normalize colors to [-1, 1] range
            vmax = max(abs(minimum(values_filtered)), abs(maximum(values_filtered)))
            vmax = max(vmax, 1e-10)  # Avoid division by zero
            
            # Create scatter plot
            GLMakie.scatter!(ax, xs_filtered, ys_filtered, zs_filtered,
                            color=values_filtered,
                            colormap=colormap,
                            colorrange=(-vmax, vmax),
                            markersize=slider_markersize.value[],
                            transparency=false)
            
            pct = 100 * n_points / length(values)
            println("Displayed $n_points points ($(round(pct, digits=1))% of total)")
        finally
            updating[] = false
        end
    end
    
    # Update on any change
    GLMakie.on(current_snapshot_idx) do _
        update_scatter!()
    end
    
    GLMakie.on(current_moment) do _
        update_scatter!()
    end
    
    GLMakie.on(slider_threshold.value) do _
        update_scatter!()
    end
    
    GLMakie.on(slider_markersize.value) do _
        update_scatter!()
    end
    
    GLMakie.on(slider_subsample.value) do _
        update_scatter!()
    end
    
    GLMakie.on(toggle_positive.active) do _
        update_scatter!()
    end
    
    GLMakie.on(toggle_negative.active) do _
        update_scatter!()
    end
    
    # Colorbar
    cb = GLMakie.Colorbar(fig[1:3, 4],
                         label="Moment Value",
                         colormap=colormap,
                         limits=(-1.0, 1.0),
                         width=25)
    
    # Statistics panel
    stats_text = GLMakie.@lift(begin
        snap_idx = $current_snapshot_idx
        S_field = snapshots[snap_idx].S
        moment_data = get_standardized_moment(S_field, $current_moment)
        name = $current_moment
        thresh = slider_threshold.value[]
        sub = slider_subsample.value[]
        
        n_above_thresh = sum(abs.(moment_data) .> thresh)
        pct_above = 100 * n_above_thresh / length(moment_data)
        
        @sprintf("""
        %s Statistics:
          Min:    %+.4f
          Max:    %+.4f
          Mean:   %+.4f
          Std:    %.4f
          
        Filtering:
          |%s| > %.2f: %d points (%.1f%%)
          Subsample: every %d points
          
        Time:
          t = %.4f
          Snapshot %d/%d
        """, 
        name,
        minimum(moment_data), maximum(moment_data),
        sum(moment_data)/length(moment_data),
        sqrt(sum((moment_data .- sum(moment_data)/length(moment_data)).^2)/length(moment_data)),
        name, thresh, n_above_thresh, pct_above,
        sub,
        snapshots[snap_idx].t,
        snap_idx, length(snapshots))
    end)
    
    controls[17, 1] = GLMakie.Label(fig, stats_text, fontsize=10, 
                                    halign=:left, tellwidth=false)
    
    # Initial plot
    update_scatter!()
    
    println("\n" * "="^70)
    println("SCATTERPLOT VIEWER READY!")
    println("="^70)
    println("\nControls:")
    println("  • Use time slider to step through snapshots")
    println("  • Click moment buttons to switch quantities")
    println("  • Adjust threshold to filter small values")
    println("  • Adjust point size for visibility")
    println("  • Subsample for better performance")
    println("  • Toggle positive/negative to filter by sign")
    println("  • Mouse: drag to rotate, scroll to zoom")
    println("\nPhysical Interpretation:")
    println("  S110: u-v velocity correlation (xy shear)")
    println("  S101: u-w velocity correlation (xz shear)")
    println("  S011: v-w velocity correlation (yz shear)")
    println("  S220, S202, S022: Temperature/energy anisotropy")
    println("  S300, S030, S003: Directional skewness")
    println("  S400, S040, S004: Directional kurtosis (Gaussian=3)")
    println("\nPress Enter to close.")
    println("="^70)
    
    display(fig)
    readline()
    
    return fig
end

export interactive_standardized_scatter

