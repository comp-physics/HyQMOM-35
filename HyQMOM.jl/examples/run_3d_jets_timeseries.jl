"""
3D Crossing Jets with Time-Series Visualization

Demonstrates:
- Full parameter control (code defaults + command-line overrides)
- 3D initial conditions with cubic jet regions
- Time-series snapshot collection
- Interactive 3D visualization over time
- High-resolution PNG export from interactive viewer
- Automatic MPI support (serial or parallel)

Requirements:
  GLMakie for interactive visualization:
    julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'

Usage:
  # Serial with interactive viewer
  julia --project=. examples/run_3d_jets_timeseries.jl
  julia --project=. examples/run_3d_jets_timeseries.jl --Np 60 --tmax 0.1
  
  # MPI parallel
  mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --Np 100
  mpiexec -n 8 julia --project=. examples/run_3d_jets_timeseries.jl --Np 120 --snapshot-interval 5
  
  # For headless systems/clusters (no X11/display) - just run simulation
  julia --project=. examples/run_3d_jets_timeseries.jl --no-viz true
  mpiexec -n 4 julia --project=. examples/run_3d_jets_timeseries.jl --no-viz true
  
  # With interactive viewer - click 'Save PNG' button to export high-resolution images
"""

# Check if visualization should be disabled (for headless systems/CI)
const VIZ_DISABLED = "--no-viz" in ARGS && 
    (idx = findfirst(x -> x == "--no-viz", ARGS)) !== nothing &&
    idx < length(ARGS) &&
    lowercase(ARGS[idx + 1]) in ["true", "t", "yes", "y", "1"]

# Load GLMakie first to enable visualization support in HyQMOM
# Only load if visualization is NOT disabled
const GLMAKIE_LOADED = if VIZ_DISABLED
    false
else
    try
        import GLMakie
        true
    catch e
        @warn """
        GLMakie is not installed. Visualization will not be available.
        To install: julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'
        """
        false
    end
end


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
    Nx = 40,
    Ny = 40,
    Nz = 40,
    tmax = 0.2,
    Ma = 0.0,
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
    # With snapshots (streaming mode)
    if rank == 0
        println("\nRunning with snapshot streaming...")
    end
    
    snapshot_filename, grid = simulation_runner(params)
    
    if rank == 0 && snapshot_filename !== nothing
        println("\n" * "="^70)
        println("SIMULATION COMPLETE")
        println("="^70)
        println("Snapshots saved to: $snapshot_filename")
        println("="^70)
        
        # Launch interactive viewer (unless disabled)
        if !params.no_viz && GLMAKIE_LOADED
            println("\nLaunching Interactive Time-Series Viewer (Streaming)...")
            println("\nViewer Controls:")
            println("  • Time slider: Step through snapshots (loaded on-demand)")
            println("  • Play/Pause: Animate the time evolution")
            println("  • Quantity buttons: Switch between Density, U, V, W velocities")
            println("  • Isosurface sliders: Adjust visualization levels")
            println("  • Mouse: Rotate (drag), Zoom (scroll)")
            println("="^70)
            
            try
                using JLD2
                interactive_3d_timeseries_streaming(snapshot_filename, grid, params)
            catch e
                @warn "Viewer failed" exception=(e, catch_backtrace())
                println("\nTo view results later, use:")
                println("  julia visualize_jld2.jl $snapshot_filename")
            end
        elseif params.no_viz
            println("\n" * "="^70)
            println("VISUALIZATION DISABLED (--no-viz flag)")
            println("="^70)
            println("Skipping interactive viewer (headless mode)")
            println("To view results later:")
            println("  julia visualize_jld2.jl $snapshot_filename")
            println("="^70)
        else
            println("\n" * "="^70)
            println("GLMakie not available - skipping visualization")
            println("="^70)
            println("To view results later:")
            println("  julia visualize_jld2.jl $snapshot_filename")
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
