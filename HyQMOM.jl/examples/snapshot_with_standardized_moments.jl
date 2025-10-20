"""
Example: Running simulation with standardized moments in snapshots

This example shows how to save standardized (and/or central) moments
along with raw moments in snapshot data for interactive visualization.
"""

using MPI
using HyQMOM
using JLD2

function main()
    MPI.Init()
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    # Setup simulation parameters
    params = (
        Nx = 40,
        Ny = 40,
        Nz = 20,
        tmax = 0.1,
        Kn = 1.0,
        Ma = 1.0,
        flag2D = 0,
        CFL = 0.5,
        dx = 1.0/40,
        dy = 1.0/40,
        dz = 1.0/20,
        Nmom = 35,
        nnmax = 10000,
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
        snapshot_interval = 10,              # Save every 10 steps
        save_standardized_moments = true,   # NEW: Save S4 field
        save_central_moments = false        # NEW: Optionally save C4 field
    )
    
    if rank == 0
        println("="^60)
        println("Running simulation with standardized moment snapshots")
        println("="^60)
    end
    
    # Run simulation
    snapshots, grid = simulation_runner(params)
    
    # Save results on rank 0
    if rank == 0
        println("\nSaving results...")
        
        # Save to JLD2 file
        @save "snapshots_with_standardized.jld2" snapshots grid
        
        println("Saved $(length(snapshots)) snapshots")
        println("\nEach snapshot contains:")
        println("  - M: Raw moments ($(size(snapshots[1].M)))")
        println("  - S: Standardized moments ($(size(snapshots[1].S)))")
        println("  - t: Time")
        println("  - step: Step number")
        
        println("\nExample: Extracting S110 (velocity correlation) at final time:")
        final_snap = snapshots[end]
        S110 = get_standardized_moment(final_snap.S, "S110")
        println("  S110 range: [$(minimum(S110)), $(maximum(S110))]")
        println("  S110 mean: $(sum(S110)/length(S110))")
        
        println("\nYou can now load and visualize these snapshots:")
        println("  using JLD2")
        println("  @load \"snapshots_with_standardized.jld2\" snapshots grid")
        println("  S022 = get_standardized_moment(snapshots[5].S, \"S022\")")
        println("  # Then plot S022[:,:,10] for z-slice 10 at snapshot 5")
    end
    
    MPI.Finalize()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

