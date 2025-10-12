"""
3D Crossing Jets Simulation with Enhanced Interactive Visualization

This example runs a true 3D crossing jets simulation and launches an
interactive 3D viewer with comprehensive visualization capabilities.

Usage:
    julia --project=. examples/run_3d_with_interactive.jl [options]

Optional command-line arguments:
    Np          Grid points per dimension (default: 40)
    tmax        Maximum simulation time (default: 0.1)
    Ma          Mach number (default: 1.0)
    CFL         CFL number (default: 0.7)
    Nz          Grid points in z-direction (default: 20)
    
Example:
    julia --project=. examples/run_3d_with_interactive.jl 60 0.2 1.0 0.6 30
"""

using HyQMOM
using MPI
using Printf

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nranks = MPI.Comm_size(comm)

# Parse command-line arguments
Np = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 40
tmax = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.1
Ma = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1.0
CFL = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.7
Nz = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 20

if rank == 0
    println("\n" * "="^70)
    println("3D CROSSING JETS SIMULATION")
    println("Enhanced Interactive Visualization")
    println("="^70)
    println("\nSimulation Parameters:")
    println("  Grid:            $(Np)×$(Np)×$(Nz)")
    println("  MPI ranks:       $(nranks)")
    println("  Time:            t ∈ [0, $(tmax)]")
    println("  Mach number:     Ma = $(Ma)")
    println("  CFL:             $(CFL)")
    println("  Knudsen number:  Kn = 1.0")
    println("  3D cubic jets:   enabled")
    println("="^70)
end

# Set up simulation parameters
params = (
    Np = Np,
    Nz = Nz,
    tmax = tmax,
    Kn = 1.0,
    Ma = Ma,
    flag2D = 0,  # 3D simulation
    CFL = CFL,
    dx = 1.0/Np,
    dy = 1.0/Np,
    dz = 1.0/Nz,
    Nmom = 35,
    nnmax = 100,
    dtmax = 1e-2,
    rhol = 1.0,      # Jet density
    rhor = 0.01,     # Background density
    T = 1.0,         # Temperature
    r110 = 0.0,      # U velocity for bottom-left jet
    r101 = 0.0,      # (currently 0, but can be changed)
    r011 = 0.0,      # (currently 0, but can be changed)
    symmetry_check_interval = 100,
    homogeneous_z = false,  # TRUE 3D cubic jets
    debug_output = false,
    enable_memory_tracking = false
)

# Run simulation
if rank == 0
    println("\nRunning simulation...")
    println("Note: Edge corner correction warnings are normal for sharp boundaries")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

# Launch interactive visualization (rank 0 only)
if rank == 0 && M_final !== nothing
    println("\n" * "="^70)
    println("Simulation complete!")
    println("  Final time: $(final_time)")
    println("  Time steps: $(time_steps)")
    println("="^70)
    
    # Check if GLMakie is available
    try
        # Launch interactive 3D viewer
        interactive_3d_viewer(M_final, grid, params,
                             n_streamlines=8,
                             vector_step=4,
                             streamline_length=50,
                             iso_threshold=0.5)
    catch e
        if isa(e, UndefVarError) && e.var == :interactive_3d_viewer
            println("\n" * "="^70)
            println("WARNING: GLMakie not available")
            println("="^70)
            println("Interactive 3D visualization requires GLMakie.")
            println("To install: julia> using Pkg; Pkg.add(\"GLMakie\")")
            println("\nFalling back to static visualization...")
            println("="^70)
            
            # Fallback to static plots
            try
                xm = collect(grid.xm)
                ym = collect(grid.ym)
                zm = collect(grid.zm)
                plot_final_results(M_final, xm, ym, zm=zm, show_plots=true)
                println("\nStatic plots generated successfully.")
            catch plot_err
                @warn "Failed to generate plots" exception=(plot_err, catch_backtrace())
            end
        else
            @warn "Failed to launch interactive viewer" exception=(e, catch_backtrace())
        end
    end
end

MPI.Finalize()

if rank == 0
    println("\n" * "="^70)
    println("DONE")
    println("="^70)
end

