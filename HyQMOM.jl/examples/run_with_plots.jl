#!/usr/bin/env julia
"""
Example: Run HyQMOM simulation with visualization

This example demonstrates how to use the plotting functionality to visualize
simulation results. It runs a small simulation and generates all the final
result plots (Figures 2-12 from MATLAB).

# Usage
```bash
cd HyQMOM.jl/examples
julia --project=.. run_with_plots.jl
```

# With MPI
```bash
mpiexec -n 2 julia --project=.. run_with_plots.jl
```

# Requirements
- PyPlot must be installed
- Python matplotlib backend must be available
"""

using HyQMOM
using Printf

println("="^70)
println("HyQMOM Simulation with Visualization")
println("="^70)

# Set up parameters for a small test case with visualization
println("\nSimulation parameters:")
println("  Grid size: 40 x 40")
println("  Final time: 0.1")
println("  Knudsen number: 1.0")
println("  Mach number: 0.0")
println("  CFL number: 0.5")
println("  Visualization: Enabled")
println()

# Run simulation with plotting enabled
results = HyQMOM.run_simulation(
    Np = 40,
    tmax = 0.1,
    Kn = 1.0,
    Ma = 0.0,
    CFL = 0.5,
    verbose = true,
    enable_plots = true,      # Enable visualization
    save_figures = false      # Set to true to save figures to disk
)

# Display results summary
if !isnothing(results[:M])
    println("\nResults Summary:")
    println("  Final time: $(results[:final_time])")
    println("  Time steps: $(results[:time_steps])")
    println("  Final moment field size: $(size(results[:M]))")
    
    # Basic checks
    M = results[:M]
    println("\nPhysical Checks:")
    ρ = M[:, :, 1]
    @printf("  Min density: %.6e\n", minimum(ρ))
    @printf("  Max density: %.6e\n", maximum(ρ))
    @printf("  Mean density: %.6e\n", mean(ρ))
    
    if all(ρ .> 0)
        println("  OK All densities positive")
    else
        println("  WARNING Some densities non-positive")
    end
    
    if !any(isnan, M)
        println("  OK No NaN values in moments")
    else
        println("  WARNING NaN values detected")
    end
    
    println("\n" * "="^70)
    println("Simulation and visualization complete!")
    println("="^70)
    println("\nGenerated figures:")
    println("  Figure 2: Moment line plots along diagonal")
    println("  Figure 3: Central moment line plots")
    println("  Figure 4: Standardized moment line plots")
    println("  Figure 9: Contour plots (12 panels)")
    println("  Figure 10: C-moment contour plots (16 panels)")
    println("  Figure 11: S-moment contour plots (12 panels)")
    println("  Figure 12: Hyperbolicity plots (9 panels)")
    println("\nClose figure windows or press Ctrl+C to exit.")
    println("="^70)
    
    # Keep script running to display plots
    # User must close windows or Ctrl+C
    try
        readline()
    catch
        println("\nExiting...")
    end
else
    println("  (Results on rank 0 only)")
end

