"""
Example: Run 3D simulation with plotting enabled

This example shows how to run a 3D simulation and generate plots.
For 3D simulations, the plots show the middle z-slice (just like MATLAB).

Requirements:
- PyPlot must be installed: julia> using Pkg; Pkg.add("PyPlot")

Usage:
    julia --project=. examples/run_3d_with_plots.jl
"""

using HyQMOM

println("="^70)
println("3D HyQMOM Simulation with Plotting")
println("="^70)

# Run 3D simulation with plotting enabled
results = run_simulation(
    Np=20,              # Grid size in x and y
    Nz=4,               # Grid size in z
    tmax=0.05,          # Simulation time
    homogeneous_z=true, # Jets at all z-levels
    enable_plots=true,  # Enable visualization
    save_figures=true,  # Save figures to disk
    output_dir="plots_3d"  # Output directory
)

println("\n" * "="^70)
println("Simulation Complete!")
println("="^70)
println("Grid: $(results[:Np])×$(results[:Np])×$(results[:Nz])")
println("Final time: $(results[:final_time])")
println("Time steps: $(results[:time_steps])")
println("Moment array shape: $(size(results[:M]))")
println("="^70)

# Print z-coordinate info
if haskey(results, :zm)
    zm = results[:zm]
    println("\nZ-coordinates ($(length(zm)) points):")
    for (k, z) in enumerate(zm)
        println("  k=$k: z = $z")
    end
    
    k_mid = div(length(zm), 2) + 1
    println("\n📊 Plots show z-slice $k_mid (z = $(zm[k_mid]))")
end

println("\n✅ Figures saved to plots_3d/")
println("="^70)

