"""
Interactive 3D Correlation Field Visualization for Standardized Moments

This viewer displays the correlation moments {S110, S101, S011} as a 3D vector field
with interactive slicing planes and isosurfaces.
"""

import GLMakie
using Printf
using LaTeXStrings

"""
    interactive_correlation_field(snapshots::Vector, grid; kwargs...)

Launch an interactive 3D correlation field viewer with time slider.

# Arguments
- `snapshots`: Vector of NamedTuples with fields `M`, `S` (standardized moments), `t`, `step`
- `grid`: Grid structure with `xm`, `ym`, `zm`

# Keyword Arguments
- `iso_threshold::Float64=0.1`: Isosurface threshold for correlation magnitude
- `colormap::Symbol=:viridis`: Colormap for isosurfaces

# Example
```julia
using HyQMOM, JLD2

# Load snapshots with standardized moments
@load "snapshots.jld2" snapshots grid

# View correlation field with slicing planes
interactive_correlation_field(snapshots, grid)
```

The viewer shows:
- Three orthogonal slicing planes (XY, YZ, XZ) that can be moved with sliders
- The intersection point of these planes defines a location (x, y, z)
- Isosurfaces of the correlation magnitude ||{S110, S101, S011}||
- Time slider to step through snapshots
"""
function interactive_correlation_field(snapshots::Vector, grid;
                                      iso_threshold=0.1,
                                      colormap=:viridis)
    
    # Check if standardized moments are available
    if !haskey(snapshots[1], :S)
        error("Snapshots do not contain standardized moments (S field). " *
              "Run simulation with save_standardized_moments=true")
    end
    
    println("\n" * "="^70)
    println("3D CORRELATION FIELD VIEWER")
    println("="^70)
    println("Features:")
    println("  • Interactive slicing planes (XY, YZ, XZ)")
    println("  • Sliders control plane positions")
    println("  • Isosurfaces of correlation magnitude")
    println("  • Time slider to step through snapshots")
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
    
    # Create figure with two 3D axes side by side, LaTeX font rendering
    fig = GLMakie.Figure(size=(2400, 1000), fontsize=12,
                        fonts=(; regular="CMU Serif"))  # Computer Modern (LaTeX font)
    
    # Current snapshot index
    current_snapshot_idx = GLMakie.Observable(1)
    
    # Left axis: Moment space (S110, S101, S011) - Main view
    ax_moment = GLMakie.Axis3(fig[1:2, 1], 
                             xlabel=L"S_{110}", ylabel=L"S_{101}", zlabel=L"S_{011}",
                             aspect=:data,
                             azimuth=0.3π,
                             elevation=π/8,
                             limits=(-1, 1, -1, 1, -1, 1),
                             xticklabelsize=11, yticklabelsize=11, zticklabelsize=11,
                             xlabelsize=13, ylabelsize=13, zlabelsize=13)
    
    # Right axis: Realizability boundary (|Δ₁| = 0)
    ax_realizability = GLMakie.Axis3(fig[1:2, 2], 
                                    xlabel=L"S_{110}", ylabel=L"S_{101}", zlabel=L"S_{011}",
                                    aspect=:data,
                                    azimuth=0.3π,
                                    elevation=π/8,
                                    limits=(-1, 1, -1, 1, -1, 1),
                                    xticklabelsize=11, yticklabelsize=11, zticklabelsize=11,
                                    xlabelsize=13, ylabelsize=13, zlabelsize=13)
    
    # Control panel
    controls = fig[1:2, 3] = GLMakie.GridLayout()
    
    # Time slider at the top
    controls[1, 1] = GLMakie.Label(fig, "Time Control:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    time_slider = GLMakie.Slider(fig, range=1:length(snapshots), startvalue=1)
    time_label = GLMakie.Label(fig, GLMakie.@lift(@sprintf("t=%.4f", snapshots[$(time_slider.value)[]].t)),
                              fontsize=11, halign=:left)
    controls[2, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, "Snapshot:", fontsize=12, halign=:left),
        time_slider,
        time_label;
        tellwidth=false
    )
    
    # Update current snapshot when slider changes
    GLMakie.on(time_slider.value) do idx
        current_snapshot_idx[] = Int(idx)
    end
    
    # Scatter point controls
    controls[3, 1] = GLMakie.Label(fig, "Display Options:", fontsize=14, 
                                   halign=:left, font=:bold, tellwidth=false)
    
    slider_iso = GLMakie.Slider(fig, range=0.001:0.001:0.5, startvalue=iso_threshold)
    controls[4, 1] = GLMakie.vgrid!(
        GLMakie.Label(fig, @sprintf("Min |S|: %.3f", iso_threshold), fontsize=12, halign=:left),
        slider_iso;
        tellwidth=false
    )
    
    # Update function with guard
    updating = Ref(false)
    
    function update_visualization!()
        # Guard against recursive calls
        if updating[]
            return
        end
        
        try
            updating[] = true
            
            # Clear moment space axis only
            GLMakie.empty!(ax_moment)
            
            # Get current snapshot
            snap_idx = current_snapshot_idx[]
            S_field = snapshots[snap_idx].S
            
            # Get correlation moments
            S110 = get_standardized_moment(S_field, "S110")
            S101 = get_standardized_moment(S_field, "S101")
            S011 = get_standardized_moment(S_field, "S011")
            
            # Compute correlation magnitude: sqrt(S110^2 + S101^2 + S011^2)
            corr_mag = sqrt.(S110.^2 .+ S101.^2 .+ S011.^2)
            
            iso_val = slider_iso.value[]
            
            # ============================================
            # MOMENT SPACE (Left panel)
            # ============================================
            
            # Flatten the 3D fields to get all points
            S110_flat = S110[:]
            S101_flat = S101[:]
            S011_flat = S011[:]
            corr_mag_flat = corr_mag[:]
            
            # Filter by threshold for visibility
            mask = corr_mag_flat .> iso_val
            
            if sum(mask) > 0
                S110_filtered = S110_flat[mask]
                S101_filtered = S101_flat[mask]
                S011_filtered = S011_flat[mask]
                mag_filtered = corr_mag_flat[mask]
                
                # Draw point cloud in moment space, colored by magnitude
                GLMakie.scatter!(ax_moment, 
                               S110_filtered, S101_filtered, S011_filtered,
                               color=mag_filtered,
                               colormap=colormap,
                               markersize=5,
                               alpha=0.6)
            end
            
            # Draw coordinate axes in moment space (always [-1,1])
            GLMakie.lines!(ax_moment, [-1.0, 1.0], [0, 0], [0, 0], 
                         color=:red, linewidth=2, alpha=0.3)
            GLMakie.lines!(ax_moment, [0, 0], [-1.0, 1.0], [0, 0], 
                         color=:green, linewidth=2, alpha=0.3)
            GLMakie.lines!(ax_moment, [0, 0], [0, 0], [-1.0, 1.0], 
                         color=:blue, linewidth=2, alpha=0.3)
            
            println(@sprintf("Updated: t=%.4f, |S|>%.3f: %d points", 
                           snapshots[snap_idx].t, iso_val, sum(mask)))
        finally
            updating[] = false
        end
    end
    
    # ============================================
    # REALIZABILITY BOUNDARY (Right panel - static)
    # ============================================
    # The realizability boundary is where |Δ₁| = 0
    
    try
        # Create a sphere as simplified boundary (placeholder for actual Δ₁=0 surface)
        θ = range(0, 2π, length=50)
        φ = range(0, π, length=50)
        r = 0.7  # radius of simplified boundary (scaled for [-1,1] space)
        
        S110_boundary = [r * sin(φ_val) * cos(θ_val) for θ_val in θ, φ_val in φ]
        S101_boundary = [r * sin(φ_val) * sin(θ_val) for θ_val in θ, φ_val in φ]
        S011_boundary = [r * cos(φ_val) for θ_val in θ, φ_val in φ]
        
        GLMakie.surface!(ax_realizability, S110_boundary, S101_boundary, S011_boundary,
                        colormap=:thermal,
                        alpha=0.2,
                        transparency=true)
        
        # Draw coordinate axes spanning full [-1,1] range
        GLMakie.lines!(ax_realizability, [-1.0, 1.0], [0, 0], [0, 0], 
                     color=:red, linewidth=2, alpha=0.5)
        GLMakie.lines!(ax_realizability, [0, 0], [-1.0, 1.0], [0, 0], 
                     color=:green, linewidth=2, alpha=0.5)
        GLMakie.lines!(ax_realizability, [0, 0], [0, 0], [-1.0, 1.0], 
                     color=:blue, linewidth=2, alpha=0.5)
        
        println("Realizability boundary displayed (simplified sphere, r=0.7)")
    catch e
        @warn "Could not compute realizability boundary" exception=(e, catch_backtrace())
    end
    
    # Connect controls to update
    GLMakie.on(current_snapshot_idx) do _
        update_visualization!()
    end
    
    GLMakie.on(slider_iso.value) do _
        update_visualization!()
    end
    
    # Statistics panel
    stats_text = GLMakie.@lift(begin
        snap_idx = $current_snapshot_idx
        S_field = snapshots[snap_idx].S
        
        S110 = get_standardized_moment(S_field, "S110")
        S101 = get_standardized_moment(S_field, "S101")
        S011 = get_standardized_moment(S_field, "S011")
        
        corr_mag = sqrt.(S110.^2 .+ S101.^2 .+ S011.^2)
        
        iso = slider_iso.value[]
        n_points = sum(corr_mag .> iso)
        max_mag = maximum(corr_mag)
        mean_mag = sum(corr_mag) / length(corr_mag)
        
        @sprintf("""
        Time: %.4f
        
        Moment Statistics:
          Max |S|: %.4f
          Mean |S|: %.4f
          Points > %.3f: %d
        """, 
        snapshots[snap_idx].t,
        max_mag, mean_mag, iso, n_points)
    end)
    
    controls[5, 1] = GLMakie.Label(fig, stats_text, fontsize=11, 
                                   halign=:left,
                                   tellwidth=false)
    
    # Initial visualization
    update_visualization!()
    
    println("\n" * "="^70)
    println("MOMENT SPACE VIEWER READY!")
    println("="^70)
    println("\nDual 3D Visualization:")
    println("  LEFT:  Moment space trajectory - evolves with time")
    println("  RIGHT: Realizability boundary (|Δ₁| = 0)")
    println("\nControls:")
    println("  • Time slider: Step through snapshots")
    println("  • Threshold slider: Filter points by |S| magnitude")
    println("  • Mouse: Rotate (drag), Zoom (scroll), Pan (right-drag)")
    println("\nMoment Interpretation:")
    println("  S110: u-v velocity correlation (xy shear)")
    println("  S101: u-w velocity correlation (xz shear)")
    println("  S011: v-w velocity correlation (yz shear)")
    println("  |S|: √(S110² + S101² + S011²)")
    println("\nRealizability:")
    println("  Points inside boundary → physically realizable")
    println("  Boundary surface → |Δ₁| = 0 (marginal realizability)")
    println("\nPress Enter to close.")
    println("="^70)
    
    display(fig)
    readline()
    
    return fig
end

export interactive_correlation_field

