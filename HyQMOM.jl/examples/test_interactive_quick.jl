"""
Quick test of interactive 3D viewer (very short simulation)
"""

using HyQMOM
using MPI

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Very short simulation for testing
params = (
    Np = 20,
    Nz = 10,
    tmax = 0.01,  # Very short
    Kn = 1.0,
    Ma = 0.5,
    flag2D = 0,
    CFL = 0.5,
    dx = 1.0/20,
    dy = 1.0/20,
    dz = 1.0/10,
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
    println("Running quick test simulation (Np=20, Nz=10, tmax=0.01)...")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

if rank == 0 && M_final !== nothing
    println("\nSimulation complete. Launching interactive viewer...")
    
    try
        interactive_3d_viewer(M_final, grid, params,
                             n_streamlines=5,
                             vector_step=3,
                             streamline_length=30,
                             iso_threshold=0.5)
    catch e
        println("\nERROR launching viewer:")
        println(e)
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
end

MPI.Finalize()

