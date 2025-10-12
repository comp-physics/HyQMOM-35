"""
Example: Run 3D crossing jets simulation with comprehensive visualization

This example demonstrates a proper 3D crossing jets simulation with:
- Cubic jet regions (not just extruded squares)
- Jets move in x-y plane (W=0 by default)
- Comprehensive 3D visualization including:
  * Multi-slice contour plots
  * 3D isosurfaces
  * Centerline profiles showing z-variation
  * Realizability diagnostics

The simulation creates two cubic jets:
- Bottom jet: Moving in +x, +y direction (diagonal)
- Top jet: Moving in -x, -y direction (opposite diagonal)
- Both jets confined to cubic regions centered at z=0

Requirements:
- PyPlot must be installed: julia> using Pkg; Pkg.add("PyPlot")
- mplot3d toolkit (part of matplotlib, accessed via PyCall)

Usage:
    julia --project=. examples/run_3d_crossing_jets.jl
    
    # Or with MPI (multiple ranks):
    mpiexecjl -n 4 julia --project=. examples/run_3d_crossing_jets.jl
"""

using HyQMOM
using MPI

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    println("="^70)
    println("3D Crossing Jets Simulation")
    println("="^70)
end

# Simulation parameters
params = (
    Np = 40,              # Grid resolution in x-y
    Nz = 20,              # Grid resolution in z
    tmax = 0.02,           # Simulation time
    Kn = 1.0,             # Knudsen number
    Ma = 1.0,             # Mach number (Ma=1 is considered an easy case)
    flag2D = 0,           # 3D simulation
    CFL = 0.7,            # CFL number
    dx = 1.0/40,          # Grid spacing (matches Np)
    dy = 1.0/40,
    dz = 1.0/20,
    Nmom = 35,            # Number of moments
    nnmax = 100000,       # Maximum number of time steps
    dtmax = 1e-2,         # Maximum time step
    
    # Initial condition parameters
    rhol = 1.0,           # Jet density
    rhor = 0.01,          # Background density
    T = 1.0,              # Temperature
    r110 = 0.0,           # Correlation coefficients
    r101 = 0.0,
    r011 = 0.0,
    
    # Diagnostic parameters
    symmetry_check_interval = 100,
    homogeneous_z = false,  
    debug_output = false,
    enable_memory_tracking = false
)

# Run simulation
if rank == 0
    println("\nRunning 3D crossing jets simulation...")
    println("  Grid: $(params.Np)×$(params.Np)×$(params.Nz)")
    println("  Jet geometry: Cubic regions (10% of domain in each direction)")
    println("  Jet velocities: U,V diagonal in x-y plane, W=0")
    println("  tmax: $(params.tmax), Ma: $(params.Ma), Kn: $(params.Kn)")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

# Visualization (only on rank 0)
if rank == 0
    println("\n" * "="^70)
    println("Simulation Complete!")
    println("="^70)
    println("Final time: $(final_time)")
    println("Time steps: $(time_steps)")
    
    if M_final !== nothing
        println("Result array shape: $(size(M_final))")
        
        # Create output directory
        output_dir = "plots_3d_crossing_jets"
        mkpath(output_dir)
        
        # Try interactive 3D viewer first (if GLMakie available)
        interactive_available = try
            isdefined(HyQMOM, :interactive_3d_viewer)
        catch
            false
        end
        
        if interactive_available
            println("\n✨ Launching Enhanced Interactive 3D Viewer (DEFAULT)...")
            println("\nInteractive features:")
            println("  • 18 quantities: Density, Speed, Pressure, Temperature,")
            println("                   U/V/W, C-moments, Deltas, H-moments")
            println("  • Adjustable slice planes, vectors, streamlines, isosurfaces")
            println("  • Mouse: drag to rotate, scroll to zoom")
            
            try
                interactive_3d_viewer(M_final, grid, params,
                                     n_streamlines=8,
                                     vector_step=4,
                                     streamline_length=50,
                                     iso_threshold=0.5)
            catch e
                @warn "Interactive viewer failed, falling back to static plots" exception=(e, catch_backtrace())
                interactive_available = false
            end
        end
        
        # Fallback to static plots if needed
        if !interactive_available
            println("\n📊 Generating static plots...")
            println("(Install GLMakie for interactive: Pkg.add(\"GLMakie\"))")
            println("Output directory: $(output_dir)/")
            
            try
                plot_final_results(
                    M_final,
                    grid.xm,
                    grid.ym,
                    params.Np,
                    params.Nmom,
                    save_figures=true,
                    output_dir=output_dir,
                    zm=grid.zm,
                    Nz=params.Nz
                )
                
                println("\n✅ All figures saved to $(output_dir)/")
                println("\nGenerated figures:")
                println("  Standard 2D plots (middle z-slice):")
                println("    - fig01_moments.png through fig07_hyperbolicity.png")
                println("\n  3D-specific plots:")
                println("    - fig101-110: slices, profiles, isosurface")
                
            catch e
                @error "Visualization failed" exception=e
            end
        end
    end
    
    println("\n" * "="^70)
    println("Analysis Tips:")
    println("="^70)
    println("1. Check density_slices to verify cubic jet structure")
    println("2. Verify U,V velocities show diagonal motion in x-y")
    println("3. Confirm W velocity is near zero (2D motion in 3D space)")
    println("4. Examine centerline_profiles for z-variation")
    println("5. Check realizability metrics (Δ₁>0, H₂₀₀≥0)")
    println("6. 3D isosurface should show two distinct cubic jets")
    println("="^70)
end

# Finalize MPI
MPI.Finalize()

