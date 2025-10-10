"""
Main entry point for HyQMOM 3D MPI-parallel simulation.

This script sets up parameters and runs the crossing jets problem.

# Usage
```bash
# Single rank
julia src/main.jl

# MPI parallel
mpiexec -n 4 julia --project src/main.jl

# With custom parameters
mpiexec -n 8 julia --project src/main.jl --Np 240 --tmax 0.04
```

# Command-line Arguments
- `--Np`: Grid size (default: 120)
- `--tmax`: Maximum simulation time (default: 0.02)
- `--Ma`: Mach number (default: 0.0, matching MATLAB)
- `--Kn`: Knudsen number (default: 1.0, matching MATLAB)
- `--CFL`: CFL number (default: 0.5, matching MATLAB)
- `--output`: Output file name (default: "results.jld2")
"""

using MPI
using Printf
using JLD2
using HyQMOM

function parse_args()
    # Default parameters (matching MATLAB test case defaults)
    params = Dict{Symbol,Any}(
        :Np => 120,
        :tmax => 0.02,
        :Ma => 0.0,      # MATLAB default (was 2.0)
        :Kn => 1.0,      # MATLAB default (was 0.01)
        :CFL => 0.5,     # MATLAB default (was 0.9)
        :flag2D => 0,
        :output => "results.jld2",
        :enable_plots => true,
        :save_figures => false,
        :output_dir => "."
    )
    
    # Parse command-line arguments
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--Np" && i < length(ARGS)
            params[:Np] = parse(Int, ARGS[i+1])
            i += 2
        elseif arg == "--tmax" && i < length(ARGS)
            params[:tmax] = parse(Float64, ARGS[i+1])
            i += 2
        elseif arg == "--Ma" && i < length(ARGS)
            params[:Ma] = parse(Float64, ARGS[i+1])
            i += 2
        elseif arg == "--Kn" && i < length(ARGS)
            params[:Kn] = parse(Float64, ARGS[i+1])
            i += 2
        elseif arg == "--CFL" && i < length(ARGS)
            params[:CFL] = parse(Float64, ARGS[i+1])
            i += 2
        elseif arg == "--output" && i < length(ARGS)
            params[:output] = ARGS[i+1]
            i += 2
        elseif arg == "--2D"
            params[:flag2D] = 1
            i += 1
        elseif arg == "--enable-plots"
            params[:enable_plots] = true
            i += 1
        elseif arg == "--no-plots"
            params[:enable_plots] = false
            i += 1
        elseif arg == "--save-figures"
            params[:save_figures] = true
            i += 1
        elseif arg == "--output-dir" && i < length(ARGS)
            params[:output_dir] = ARGS[i+1]
            i += 2
        elseif arg == "--help" || arg == "-h"
            println("""
            Usage: julia src/main.jl [OPTIONS]
            
            Options:
              --Np N             Grid size (default: 120)
              --tmax T           Maximum simulation time (default: 0.02)
              --Ma M             Mach number (default: 0.0)
              --Kn K             Knudsen number (default: 1.0)
              --CFL C            CFL number (default: 0.5)
              --2D               Run 2D simulation (default: 3D)
              --output FILE      Output file name (default: results.jld2)
              --enable-plots     Enable visualization (default: true)
              --no-plots         Disable visualization
              --save-figures     Save figures to disk (default: false)
              --output-dir DIR   Directory for saved figures (default: .)
              --help, -h         Show this help message
            
            Examples:
              # Run with default visualization
              julia --project main.jl --Np 40 --tmax 0.1
              
              # Run without visualization
              julia --project main.jl --Np 40 --tmax 0.1 --no-plots
              
              # Run and save figures
              julia --project main.jl --Np 40 --tmax 0.1 --save-figures --output-dir ./figures
              
              # MPI parallel
              mpiexec -n 4 julia --project main.jl --Np 240 --tmax 0.04
            """)
            exit(0)
        else
            i += 1
        end
    end
    
    return params
end

function main()
    # Initialize MPI
    MPI.Init()
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    # Parse command-line arguments
    args = parse_args()
    
    # Setup simulation parameters
    Np = args[:Np]
    tmax = args[:tmax]
    Ma = args[:Ma]
    Kn = args[:Kn]
    CFL = args[:CFL]
    flag2D = args[:flag2D]
    output_file = args[:output]
    enable_plots = args[:enable_plots]
    save_figures = args[:save_figures]
    output_dir = args[:output_dir]
    
    # Derived parameters
    dx = 1.0 / Np
    dy = 1.0 / Np
    Nmom = 35
    nnmax = 100000
    dtmax = 0.001
    
    # Initial condition parameters (crossing jets, matching MATLAB)
    rhol = 1.0    # High density in jets
    rhor = 0.01   # Low density in background (MATLAB default, was 0.5)
    T = 1.0       # Temperature
    r110 = 0.0    # Correlation coefficients
    r101 = 0.0
    r011 = 0.0
    
    # Diagnostic parameters
    symmetry_check_interval = 10
    
    # Package parameters into named tuple
    params = (Np=Np, tmax=tmax, Kn=Kn, Ma=Ma, flag2D=flag2D, CFL=CFL,
              dx=dx, dy=dy, Nmom=Nmom, nnmax=nnmax, dtmax=dtmax,
              rhol=rhol, rhor=rhor, T=T, r110=r110, r101=r101, r011=r011,
              symmetry_check_interval=symmetry_check_interval,
              enable_memory_tracking=false,
              debug_output=false)
    
    if rank == 0
        println("="^70)
        println("HyQMOM 3D MPI-Parallel Simulation")
        println("="^70)
        println("Configuration:")
        println("  Grid size: $(Np)x$(Np)")
        println("  MPI ranks: $(nprocs)")
        println("  Max time: $(tmax)")
        println("  Mach number: $(Ma)")
        println("  Knudsen number: $(Kn)")
        println("  CFL number: $(CFL)")
        println("  2D mode: $(flag2D == 1 ? "Yes" : "No")")
        println("  Output file: $(output_file)")
        println("  Visualization: $(enable_plots ? "Enabled" : "Disabled")")
        if enable_plots && save_figures
            println("  Save figures: Yes (to $output_dir)")
        end
        println("="^70)
        flush(stdout)
    end
    
    # Run simulation
    start_time = time()
    M_final, final_time, time_steps, grid_out = simulation_runner(params)
    total_time = time() - start_time
    
    # Generate plots if requested (only on rank 0)
    if enable_plots && rank == 0 && !isnothing(M_final)
        println("="^70)
        println("Generating visualization...")
        println("  enable_plots = $enable_plots")
        println("  save_figures = $save_figures")
        println("  output_dir = $output_dir")
        println("="^70)
        plot_final_results(M_final, grid_out.xm, grid_out.ym, Np, 35;
                         save_figures=save_figures, output_dir=output_dir)
    elseif !enable_plots && rank == 0
        println("Visualization disabled (--no-plots)")
    elseif rank != 0
        # Other ranks don't plot
    elseif isnothing(M_final)
        @warn "M_final is nothing - cannot generate plots"
    end
    
    # Save results (only rank 0)
    if rank == 0
        println("="^70)
        println("Simulation complete!")
        println("  Final time: $(final_time)")
        println("  Time steps: $(time_steps)")
        println("  Wall time: $(total_time) seconds")
        println("  Time per step: $(total_time/time_steps) seconds")
        println("="^70)
        
        # Save results to file
        println("Saving results to $(output_file)...")
        jldsave(output_file;
                M_final=M_final,
                final_time=final_time,
                time_steps=time_steps,
                grid=grid_out,
                params=params)
        println("Results saved successfully!")
    end
    
    # Finalize MPI
    MPI.Finalize()
    
    return 0
end

# Run main if this is the entry script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
