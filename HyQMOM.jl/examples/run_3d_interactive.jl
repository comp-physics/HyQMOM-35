"""
3D Crossing Jets with Interactive Visualization

This launches an interactive 3D viewer using GLMakie that allows you to:
- Rotate, zoom, and pan the 3D view
- Toggle between different quantities (density, U, V, W)
- Adjust slice plane positions interactively
- Change transparency and isovalues in real-time
- View volume rendering and isosurfaces

Usage:
    julia --project=. examples/run_3d_interactive.jl
"""

using HyQMOM
using MPI
using GLMakie

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 0
    println("="^70)
    println("3D Crossing Jets - Interactive Visualization")
    println("="^70)
end

# Run simulation (same parameters as before)
params = (
    Np = 40,
    Nz = 20,
    tmax = 0.1,
    Kn = 1.0,
    Ma = 1.0,
    flag2D = 0,
    CFL = 0.7,
    dx = 1.0/40,
    dy = 1.0/40,
    dz = 1.0/20,
    Nmom = 35,
    nnmax = 100000,
    dtmax = 1e-2,
    rhol = 1.0,
    rhor = 0.01,
    T = 1.0,
    r110 = 0.0,
    r101 = 0.0,
    r011 = 0.0,
    symmetry_check_interval = 100,
    homogeneous_z = false,  # TRUE 3D cubic jets
    debug_output = false,
    enable_memory_tracking = false
)

if rank == 0
    println("\nRunning simulation...")
    println("  Grid: $(params.Np)×$(params.Np)×$(params.Nz)")
    println("  Ma: $(params.Ma), homogeneous_z: $(params.homogeneous_z)")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

# Launch interactive viewer (only on rank 0)
if rank == 0 && M_final !== nothing
    println("\n" * "="^70)
    println("Simulation Complete!")
    println("  Final time: $(final_time)")
    println("  Steps: $(time_steps)")
    println("="^70)
    
    println("\nLaunching interactive 3D viewer...")
    println("Close the window when done exploring.")
    println("="^70)
    
    # Include the interactive visualization code
    include("../src/visualization/interactive_3d.jl")
    
    # Launch viewer
    fig = interactive_3d_viewer(
        M_final,
        grid.xm,
        grid.ym,
        grid.zm,
        params.Np,
        params.Nz,
        params.Nmom,
        quantities=["density", "U", "V", "W"]
    )
    
    # Keep window open
    println("\nInteractive viewer is running...")
    println("Press Ctrl+C when done.")
    wait(fig.scene)
end

MPI.Finalize()

