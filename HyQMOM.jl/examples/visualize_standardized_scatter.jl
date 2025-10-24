"""
Example: Interactive 3D Scatterplot Visualization of Standardized Moments

This example shows how to visualize standardized moments as 3D point clouds
using GLMakie's interactive scatterplot viewer.
"""

using MPI
using HyQMOM
using JLD2

function main()
    # Check if GLMakie is available
    try
        import GLMakie
    catch
        println("ERROR: GLMakie not installed. Install with:")
        println("  using Pkg; Pkg.add(\"GLMakie\")")
        return
    end
    
    # Initialize MPI
    MPI.Init()
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    # Run a short simulation with standardized moments
    if rank == 0
        println("="^70)
        println("STANDARDIZED MOMENTS SCATTERPLOT EXAMPLE")
        println("="^70)
        println("Running short simulation with snapshot saving...")
    end
    
    # Smaller grid for quick demo
    params = (
        Nx = 30,
        Ny = 30,
        Nz = 15,
        tmax = 0.05,  # Short simulation
        Kn = 1.0,
        Ma = 1.0,
        flag2D = 0,
        CFL = 0.5,
        dx = 1.0/30,
        dy = 1.0/30,
        dz = 1.0/15,
        Nmom = 35,
        nnmax = 20000000,
        dtmax = 0.01,
        rhol = 1.0,
        rhor = 0.01,
        T = 1.0,
        r110 = 0.0,
        r101 = 0.0,
        r011 = 0.0,
        symmetry_check_interval = 10,
        homogeneous_z = true,
        debug_output = false,
        snapshot_interval = 5,              # Save every 5 steps
        save_standardized_moments = true   # IMPORTANT: Enable S field
    )
    
    # Run simulation
    snapshots, grid = simulation_runner(params)
    
    MPI.Finalize()
    
    # Visualize on rank 0
    if rank == 0
        println("\n" * "="^70)
        println("Simulation complete! Launching interactive viewer...")
        println("="^70)
        
        # Choose a snapshot to visualize (last one)
        snap = snapshots[end]
        
        println("\nViewing snapshot $(length(snapshots)): t=$(snap.t), step=$(snap.step)")
        println("Available fields: $(keys(snap))")
        
        # Launch interactive scatterplot viewer
        println("\nLaunching interactive 3D scatterplot viewer...")
        println("You can:")
        println("  • Click buttons to switch between moments (S110, S022, etc.)")
        println("  • Adjust threshold slider to filter small values")
        println("  • Change point size for better visibility")
        println("  • Subsample to improve performance")
        println("  • Toggle positive/negative values")
        println("  • Rotate/zoom with mouse")
        
        interactive_standardized_scatter(snap, grid;
                                        threshold=0.15,
                                        subsample=1,
                                        markersize=4.0,
                                        colormap=:RdBu)
        
        println("\nDone!")
    end
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

