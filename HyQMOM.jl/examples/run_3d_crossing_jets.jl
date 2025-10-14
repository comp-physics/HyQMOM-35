"""
3D Crossing Jets with Static Visualization (PyPlot)

Demonstrates 3D simulation with static matplotlib/PyPlot visualization.
Use this if you prefer static plots or don't have GLMakie installed.

For interactive 3D visualization, use: run_3d_jets_timeseries.jl

Usage:
  # Serial
  julia --project=. examples/run_3d_crossing_jets.jl
  julia --project=. examples/run_3d_crossing_jets.jl --Np 60 --tmax 0.1
  
  # MPI parallel
  mpiexec -n 4 julia --project=. examples/run_3d_crossing_jets.jl --Np 100

Requirements:
  - PyPlot: julia> using Pkg; Pkg.add("PyPlot")
"""

using HyQMOM
using MPI

# Load parameter parsing utilities
include("parse_params.jl")

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Parse parameters
params = parse_simulation_params(
    Np = 40,
    Nz = 20,
    tmax = 0.02,
    Ma = 0.0,
    Kn = 1.0,
    CFL = 0.7,
    homogeneous_z = false,
    snapshot_interval = 0  # No snapshots for static plotting
)

if rank == 0
    println("="^70)
    println("3D CROSSING JETS - STATIC VISUALIZATION")
    println("="^70)
end

print_params_summary(params, rank=rank, comm=comm)

# Run simulation
if rank == 0
    println("\nRunning simulation...")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

# Plot on rank 0
if rank == 0
    println("\n" * "="^70)
    println("SIMULATION COMPLETE")
    println("="^70)
    println("Final time: $final_time")
    println("Time steps: $time_steps")
    println("="^70)
    
    println("\nGenerating static plots...")
    try
        # Generate matplotlib/PyPlot visualizations
        plot_final_results(
            M_final,
            grid.xm, grid.ym, grid.zm,
            final_time
        )
        println("✓ Plots generated successfully")
    catch e
        @warn "Plotting failed" exception=(e, catch_backtrace())
        println("Note: Install PyPlot if not available: using Pkg; Pkg.add(\"PyPlot\")")
    end
end

MPI.Finalize()

if rank == 0
    println("\nDone!")
end
