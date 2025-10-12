"""
Test the improved interactive 3D viewer with working quantity selection
"""

using HyQMOM
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Quick simulation for testing
params = (
    Np = 30,
    Nz = 15,
    tmax = 0.02,
    Kn = 1.0,
    Ma = 0.7,
    flag2D = 0,
    CFL = 0.6,
    dx = 1.0/30,
    dy = 1.0/30,
    dz = 1.0/15,
    Nmom = 35,
    nnmax = 100,
    dtmax = 1e-2,
    rhol = 1.0,
    rhor = 0.01,
    T = 1.0,
    r110 = 0.0,
    r101 = 0.0,
    r011 = 0.0,
    symmetry_check_interval = 100,
    homogeneous_z = false,
    debug_output = false,
    enable_memory_tracking = false
)

if rank == 0
    println("Running test simulation (Np=30, Nz=15, tmax=0.02)...")
    println("This should take about 5-10 seconds...")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

if rank == 0 && M_final !== nothing
    println("\nSimulation complete!")
    println("Launching IMPROVED interactive viewer...")
    println("\nKey improvements:")
    println("  ✓ Quantity selection with buttons (not dropdown)")
    println("  ✓ All quantities should work: Density, U, V, W, etc.")
    println("  ✓ Better volume rendering")
    println("\nTry clicking the buttons to switch between quantities!")
    
    try
        interactive_3d_viewer_improved(M_final, grid, params,
                                      n_streamlines=6,
                                      vector_step=4,
                                      streamline_length=40,
                                      iso_threshold=0.5)
    catch e
        println("\nERROR launching improved viewer:")
        println(e)
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
end

MPI.Finalize()

