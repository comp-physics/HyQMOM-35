"""
Test the TRUE 3D volume/isosurface viewer

This viewer shows ACTUAL 3D surfaces using GLMakie's contour function,
not just dots or slice planes.
"""

using HyQMOM
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Quick simulation
params = (
    Np = 60,
    Nz = 30,
    tmax = 0.02,
    Kn = 1.0,
    Ma = 0.7,
    flag2D = 0,
    CFL = 0.6,
    dx = 1.0/60,
    dy = 1.0/60,
    dz = 1.0/30,
    Nmom = 35,
    nnmax = 1,
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
    println("Running simulation (Np=30, Nz=15, tmax=0.02)...")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

if rank == 0 && M_final !== nothing
    println("\n" * "="^70)
    println("Simulation complete!")
    println("="^70)
    println("\nLaunching VOLUME/ISOSURFACE viewer...")
    println("\nThis viewer shows TRUE 3D SURFACES:")
    println("  • THREE colored isosurface contours (blue, green, red)")
    println("  • NOT dots - actual 3D surfaces you can see from all angles")
    println("  • Adjust sliders to change isosurface levels")
    println("  • Rotate with mouse to view from different angles")
    println("\nWhat to expect:")
    println("  • Blue surface   = low density regions (30% of max)")
    println("  • Green surface  = medium density (50% of max)")
    println("  • Red surface    = high density (70% of max)")
    println("  • Cyan lines     = streamlines showing flow paths")
    println("\nTry this:")
    println("  1. Rotate the view to see the 3D structure")
    println("  2. Click 'Speed' button to see speed isosurfaces")
    println("  3. Adjust 'Iso level' sliders to explore different thresholds")
    println("  4. Toggle streamlines on/off")
    
    try
        interactive_3d_volume(M_final, grid, params,
                             n_streamlines=6,
                             vector_step=4,
                             streamline_length=40,
                             iso_levels=[0.3, 0.5, 0.7],
                             enable_volume=true)
    catch e
        println("\nERROR launching volume viewer:")
        showerror(stdout, e, catch_backtrace())
        println()
    end
end

MPI.Finalize()

