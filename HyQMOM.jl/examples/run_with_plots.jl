#!/usr/bin/env julia
"""
Example: Run HyQMOM simulation with visualization

This example demonstrates the high-level run_simulation() API with visualization.
It runs a small simulation and generates comprehensive diagnostic plots.

This is the recommended way to run simulations programmatically.
For the low-level API, see run_simple.jl.

# Usage
```bash
cd HyQMOM.jl
julia --project=. examples/run_with_plots.jl
```

# With MPI
```bash
cd HyQMOM.jl
mpiexec -n 2 julia --project=. examples/run_with_plots.jl
```

# Requirements
- PyPlot must be installed for visualization
- Python matplotlib backend must be available
- Set enable_plots=false if PyPlot is not available
"""

using HyQMOM
using Printf
using Statistics
using MPI

println("="^70)
println("HyQMOM Simulation with Visualization")
println("="^70)
println("\nThis example uses the high-level run_simulation() API.")
println("The API handles MPI initialization automatically.\n")

# Set up parameters for a small test case with visualization
println("Simulation parameters:")
println("  Grid size: 40 x 40")
println("  Final time: 0.1")
println("  Knudsen number: 1.0")
println("  Mach number: 0.0")
println("  CFL number: 0.5")
println("  Visualization: Enabled")
println()

# Run simulation using high-level API
# This automatically handles:
# - MPI initialization (if not already done)
# - Parameter validation
# - Grid setup
# - Simulation execution
# - Result gathering
# - Optional visualization
results = HyQMOM.run_simulation(
    # Grid and time
    Np = 40,              # Grid size: 40x40
    tmax = 0.1,           # Final time
    
    # Physical parameters
    Kn = 1.0,             # Knudsen number
    Ma = 0.0,             # Mach number
    CFL = 0.5,            # CFL number
    flag2D = 0,           # 0 = 3D, 1 = 2D
    
    # Execution options
    num_workers = MPI.Initialized() ? MPI.Comm_size(MPI.COMM_WORLD) : 1,
    verbose = true,       # Print progress
    
    # Visualization options
    enable_plots = true,  # Generate plots (requires PyPlot)
    save_figures = false, # Set to true to save figures to disk
    output_dir = ".",     # Directory for saved figures
    
    # Debug options
    debug_output = false  # Enable debug output
)

# Display results summary (only on rank 0)
rank = MPI.Initialized() ? MPI.Comm_rank(MPI.COMM_WORLD) : 0

if rank == 0 && !isnothing(results[:M])
    println("\n" * "="^70)
    println("Results Summary")
    println("="^70)
    println("  Final time: $(results[:final_time])")
    println("  Time steps: $(results[:time_steps])")
    println("  Grid size: $(results[:Np]) x $(results[:Np])")
    println("  Moment array: $(size(results[:M]))")
    
    # Extract moment field
    M = results[:M]
    
    # Physical checks
    println("\n" * "-"^70)
    println("Physical Validation")
    println("-"^70)
    
    # Check density (first moment)
    ρ = M[:, :, 1]
    @printf("  Density range: [%.6e, %.6e]\n", minimum(ρ), maximum(ρ))
    @printf("  Mean density:  %.6e\n", mean(ρ))
    @printf("  Std density:   %.6e\n", std(ρ))
    
    # Check for numerical issues
    checks_passed = true
    
    if any(isnan, M)
        println("  ✗ NaN values detected in moments")
        checks_passed = false
    else
        println("  ✓ No NaN values")
    end
    
    if any(isinf, M)
        println("  ✗ Inf values detected in moments")
        checks_passed = false
    else
        println("  ✓ No Inf values")
    end
    
    if any(ρ .<= 0)
        println("  ✗ Non-positive densities detected")
        checks_passed = false
    else
        println("  ✓ All densities positive")
    end
    
    # Check energy (fifth moment)
    E = M[:, :, 5]
    if any(E .<= 0)
        println("  ✗ Non-positive energies detected")
        checks_passed = false
    else
        println("  ✓ All energies positive")
    end
    
    if checks_passed
        println("\n  ✓✓✓ All physical checks passed!")
    else
        println("\n  ⚠⚠⚠ Some physical checks failed!")
    end
    
    # Mass conservation (rough estimate)
    println("\n" * "-"^70)
    println("Conservation Check")
    println("-"^70)
    dx = results[:xm][2] - results[:xm][1]
    dy = results[:ym][2] - results[:ym][1]
    total_mass = sum(ρ) * dx * dy
    @printf("  Total mass: %.6e\n", total_mass)
    @printf("  Domain:     [%.2f, %.2f] x [%.2f, %.2f]\n", 
            results[:xm][1]-dx/2, results[:xm][end]+dx/2,
            results[:ym][1]-dy/2, results[:ym][end]+dy/2)
    
    println("\n" * "="^70)
    println("Visualization Complete!")
    println("="^70)
    
    if get(ENV, "HYQMOM_SKIP_PLOTTING", "false") != "true"
        println("\nGenerated figures:")
        println("  Figure 2:  Moment line plots along diagonal")
        println("  Figure 3:  Central moment line plots")
        println("  Figure 4:  Standardized moment line plots")
        println("  Figure 9:  Contour plots (12 panels)")
        println("  Figure 10: C-moment contour plots (16 panels)")
        println("  Figure 11: S-moment contour plots (12 panels)")
        println("  Figure 12: Hyperbolicity plots (9 panels)")
        println()
        println("Close figure windows or press Ctrl+C to exit.")
        println("="^70)
        
        # Keep script running to display plots
        # User must close windows or press Ctrl+C
        try
            println("\nPress Enter to exit...")
            readline()
        catch
            println("\nExiting...")
        end
    else
        println("\nPlotting was skipped (HYQMOM_SKIP_PLOTTING is set)")
    end
    
elseif rank == 0
    println("\n(No results returned - may be running on non-root rank)")
elseif rank != 0
    # Non-root ranks don't print anything
end

