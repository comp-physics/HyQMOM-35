"""
3D Crossing Jets with Time-Series Visualization

Demonstrates:
- Full parameter control (code defaults + command-line overrides)
- 3D initial conditions with cubic jet regions
- Time-series snapshot collection
- Interactive 3D visualization over time
- Automatic MPI support (serial or parallel)

Usage:
  # Serial
  julia --project=. examples/run_3d_jets_timeseries.jl
  julia --project=. examples/run_3d_jets_timeseries.jl --Np 60 --tmax 0.1
  
  # MPI parallel
  mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100
  mpiexec -n 8 julia --project=. examples/run_3d_jets_timeseries.jl --Np 120 --snapshot-interval 5
"""

using HyQMOM
using MPI
using Printf

# Load parameter parsing utilities
include("parse_params.jl")

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Parse parameters: code defaults + command-line overrides
params = parse_simulation_params(
    # Code defaults (can be overridden with command-line args)
    Np = 40,
    Nz = 40,
    tmax = 0.2,
    Ma = 1.0,
    Kn = 1.0,
    CFL = 0.7,
    xmin = 0.0,
    xmax = 1.0,
    ymin = 0.0,
    ymax = 1.0,
    zmin = 0.0,
    zmax = 1.0,
    snapshot_interval = 1,
    homogeneous_z = false
)

# Print parameter summary
if rank == 0
    println("="^70)
    println("3D CROSSING JETS - TIME SERIES VISUALIZATION")
    println("="^70)
end

print_params_summary(params, rank=rank, comm=comm)

if rank == 0
    println("\nJet Configuration:")
    println("  True 3D cubic jets (homogeneous_z = false)")
    println("  Bottom-left jet: velocity (+U, +V) → moving ↗")
    println("  Top-right jet:   velocity (-U, -V) → moving ↙")
    println("="^70)
end

# Run simulation
if params.snapshot_interval > 0
    # With snapshots
    if rank == 0
        println("\nRunning with snapshot collection...")
    end
    
    snapshots, grid = run_simulation_with_snapshots(params; 
                                                     snapshot_interval=params.snapshot_interval)
    
    if rank == 0 && snapshots !== nothing
        println("\n" * "="^70)
        println("SIMULATION COMPLETE")
        println("="^70)
        println("Collected $(length(snapshots)) snapshots")
        println("\nSnapshot Timeline:")
        for (i, snap) in enumerate(snapshots)
            if i <= 5 || i > length(snapshots) - 5
                @printf("  %2d: t = %.4f, step = %d\n", i, snap.t, snap.step)
            elseif i == 6
                println("  ...")
            end
        end
        println("="^70)
        
        # Launch interactive viewer
        println("\nLaunching Interactive Time-Series Viewer...")
        println("\nViewer Controls:")
        println("  • Time slider: Step through snapshots")
        println("  • Play/Pause/Reset: Animate the time evolution")
        println("  • Quantity buttons: Switch between Density, U, V, W velocities")
        println("  • Isosurface sliders: Adjust visualization levels")
        println("  • Mouse: Rotate (drag), Zoom (scroll)")
        println("="^70)
        
        try
            interactive_3d_timeseries(snapshots, grid, params)
        catch e
            @warn "Viewer failed" exception=(e, catch_backtrace())
            println("Trying fallback single-frame viewer...")
            try
                interactive_3d_volume(snapshots[end].M, grid, params)
            catch e2
                @warn "Fallback viewer also failed"
            end
        end
    end
else
    # Without snapshots
    if rank == 0
        println("\nRunning without snapshots (snapshot_interval=0)...")
    end
    
    M_final, final_time, time_steps, grid = simulation_runner(params)
    
    if rank == 0
        println("\n" * "="^70)
        println("SIMULATION COMPLETE")
        println("="^70)
        println("Final time: $final_time")
        println("Time steps: $time_steps")
        println("="^70)
    end
end

MPI.Finalize()

if rank == 0
    println("\nDone!")
end
