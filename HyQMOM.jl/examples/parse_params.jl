"""
Unified parameter parsing for all examples.

This module provides:
1. Default parameters
2. Command-line argument parsing
3. Merging of defaults + file overrides + command-line overrides

Usage:
    include("parse_params.jl")
    
    # Option 1: Use all defaults
    params = parse_simulation_params()
    
    # Option 2: Override some defaults in code
    params = parse_simulation_params(
        Nx = 60,
        Ny = 60,
        Nz = 30,
        tmax = 0.1
    )
    
    # Option 3: Command-line args override everything
    # julia example.jl --Np 60 --tmax 0.1 --Ma 1.5
"""

"""
    get_default_params()

Returns a NamedTuple with all default simulation parameters.
"""
function get_default_params()
    return (
        # Grid resolution
        Nx = 40,
        Ny = 40,
        Nz = 20,
        
        # Domain extents (default: [-0.5, 0.5]³)
        xmin = -0.5,
        xmax = 0.5,
        ymin = -0.5,
        ymax = 0.5,
        zmin = -0.5,
        zmax = 0.5,
        
        # Time parameters
        tmax = 0.05,
        nnmax = 100,
        dtmax = 1e-2,
        CFL = 0.7,
        
        # Physical parameters
        Ma = 0.0,
        Kn = 1.0,
        T = 1.0,
        
        # Initial condition
        rhol = 1.0,
        rhor = 0.01,
        r110 = 0.0,
        r101 = 0.0,
        r011 = 0.0,
        
        # Numerical parameters
        Nmom = 35,
        flag2D = 0,
        
        # Diagnostics
        symmetry_check_interval = 100,
        homogeneous_z = false,
        debug_output = false,
        enable_memory_tracking = false,
        
        # Snapshot control
        snapshot_interval = 0,  # 0 = disabled
        save_standardized_moments = false,  # Save S4 field with snapshots
        save_central_moments = false,  # Save C4 field with snapshots
        
        # Initial condition configuration
        config = "crossing",  # Configuration name for custom ICs
    )
end

"""
    parse_command_line_args(args=ARGS)

Parse command-line arguments and return a dictionary of parameter overrides.

Supported syntax:
  --Np 60
  --tmax 0.1
  --Ma 1.5
  --xmin 0.0 --xmax 1.0
  --snapshot-interval 5
  --homogeneous-z true
  --help
"""
function parse_command_line_args(args=ARGS)
    overrides = Dict{Symbol, Any}()
    
    i = 1
    while i <= length(args)
        arg = args[i]
        
        if arg == "--help" || arg == "-h"
            print_help()
            exit(0)
        elseif startswith(arg, "--")
            # Remove leading "--" and convert hyphens to underscores
            param_name = Symbol(replace(arg[3:end], "-" => "_"))
            
            if i < length(args)
                value_str = args[i + 1]
                
                # Try to parse the value
                value = try_parse_value(value_str)
                overrides[param_name] = value
                i += 2
            else
                @warn "Parameter $arg requires a value"
                i += 1
            end
        else
            @warn "Unknown argument: $arg (use --help for usage)"
            i += 1
        end
    end
    
    return overrides
end

"""
    try_parse_value(str)

Attempt to parse a string as Int, Float64, or Bool, otherwise return as String.
"""
function try_parse_value(str)
    # Try Bool first
    if lowercase(str) in ["true", "t", "yes", "y", "1"]
        return true
    elseif lowercase(str) in ["false", "f", "no", "n", "0"]
        return false
    end
    
    # Try Int
    try
        return parse(Int, str)
    catch
    end
    
    # Try Float
    try
        return parse(Float64, str)
    catch
    end
    
    # Return as string
    return str
end

"""
    print_help()

Print usage information.
"""
function print_help()
    println("""
    Usage: julia example.jl [OPTIONS]
    
    Common Options:
      --Nx N                Grid resolution in x direction (default: 40)
      --Ny N                Grid resolution in y direction (default: 40)
      --Nz N                Grid resolution in z direction (default: 20)
      
      --xmin X              Domain extent x minimum (default: -0.5)
      --xmax X              Domain extent x maximum (default: 0.5)
      --ymin Y              Domain extent y minimum (default: -0.5)
      --ymax Y              Domain extent y maximum (default: 0.5)
      --zmin Z              Domain extent z minimum (default: -0.5)
      --zmax Z              Domain extent z maximum (default: 0.5)
      
      --tmax T              Maximum simulation time (default: 0.05)
      --nnmax N             Maximum time steps (default: 100)
      --CFL C               CFL number (default: 0.7)
      --dtmax D             Maximum time step (default: 1e-2)
      
      --Ma M                Mach number (default: 1.0)
      --Kn K                Knudsen number (default: 1.0)
      --T T                 Temperature (default: 1.0)
      
      --rhol R              Jet density (default: 1.0)
      --rhor R              Background density (default: 0.01)
      
      --snapshot-interval N           Save snapshot every N steps (default: 0=disabled)
      --save-standardized-moments BOOL Save S4 field with snapshots (default: false)
      --save-central-moments BOOL      Save C4 field with snapshots (default: false)
      
      --config NAME         Initial condition configuration (default: crossing)
                           Options: crossing, triple-jet, quad-jet, vertical-jet, spiral
      
      --homogeneous-z BOOL  Use z-homogeneous mode (default: false)
      --debug-output BOOL   Enable debug output (default: false)
      
      -h, --help            Show this help message
    
    Examples:
      julia example.jl
      julia example.jl --Nx 60 --Ny 60 --Nz 30
      julia example.jl --tmax 0.1 --Ma 1.5
      julia example.jl --xmin 0 --xmax 1 --ymin 0 --ymax 1 --zmin 0 --zmax 1
      julia example.jl --snapshot-interval 5
      mpiexec -n 4 julia example.jl --Nx 80 --Ny 80 --snapshot-interval 10
    """)
end

"""
    parse_simulation_params(; kwargs...)

Main function to parse simulation parameters.

Priority (highest to lowest):
  1. Command-line arguments
  2. Keyword arguments (code overrides)
  3. Default values

Returns a NamedTuple with all parameters.
"""
function parse_simulation_params(; kwargs...)
    # Start with defaults
    defaults = get_default_params()
    
    # Apply code overrides
    code_overrides = Dict{Symbol, Any}(kwargs)
    
    # Apply command-line overrides
    cli_overrides = parse_command_line_args()
    
    # Merge: defaults < code_overrides < cli_overrides
    params_dict = Dict{Symbol, Any}(pairs(defaults))
    merge!(params_dict, code_overrides)
    merge!(params_dict, cli_overrides)
    
    # Convert back to NamedTuple
    return NamedTuple(params_dict)
end

"""
    print_params_summary(params; rank=0, comm=nothing)

Print a formatted summary of simulation parameters.
Only prints on rank 0 by default.
If comm is provided, prints MPI information.
"""
function print_params_summary(params; rank=0, comm=nothing)
    if rank != 0
        return
    end
    
    println("="^70)
    println("SIMULATION PARAMETERS")
    println("="^70)
    
    # MPI information (if applicable)
    if comm !== nothing
        try
            # MPI is being used
            nprocs = MPI.Comm_size(comm)
            if nprocs > 1
                println("\nMPI Configuration:")
                println("  Ranks: $nprocs")
                println("  Decomposition: x-y plane divided among ranks")
                println("  Z-direction: No decomposition (all ranks have full z)")
                
                # Estimate local grid size (assumes uniform distribution)
                local_nx = div(params.Nx, nprocs) + (params.Nx % nprocs > 0 ? 1 : 0)
                println("  Estimated local grid per rank: ~$(local_nx)×$(params.Ny)×$(params.Nz) (uniform distribution estimate)")
            else
                println("\nExecution: Single rank (serial)")
            end
        catch
            # MPI not properly initialized or error, skip MPI info
        end
    else
        println("\nExecution: Serial (no MPI)")
    end
    
    println("\nGrid & Domain:")
    println("  Resolution: $(params.Nx)×$(params.Ny)×$(params.Nz)")
    println("  Domain: [$(params.xmin), $(params.xmax)] × [$(params.ymin), $(params.ymax)] × [$(params.zmin), $(params.zmax)]")
    
    dx = (params.xmax - params.xmin) / params.Nx
    dy = (params.ymax - params.ymin) / params.Ny
    dz = (params.zmax - params.zmin) / params.Nz
    println("  Grid spacing: dx=$(round(dx, digits=4)), dy=$(round(dy, digits=4)), dz=$(round(dz, digits=4))")
    
    println("\nTime Integration:")
    println("  tmax: $(params.tmax)")
    println("  nnmax: $(params.nnmax)")
    println("  CFL: $(params.CFL)")
    println("  dtmax: $(params.dtmax)")
    
    println("\nPhysics:")
    println("  Ma: $(params.Ma)")
    println("  Kn: $(params.Kn)")
    println("  T: $(params.T)")
    println("  Jet density: $(params.rhol), Background: $(params.rhor)")
    
    println("\nNumerics:")
    println("  Nmom: $(params.Nmom)")
    println("  flag2D: $(params.flag2D) (0=3D, 1=2D)")
    println("  homogeneous_z: $(params.homogeneous_z)")
    
    if params.snapshot_interval > 0
        println("\nSnapshots:")
        println("  Interval: every $(params.snapshot_interval) step(s)")
        if params.save_standardized_moments
            println("  Saving: Raw moments (M) + Standardized moments (S)")
        end
        if params.save_central_moments
            println("  Saving: Raw moments (M) + Central moments (C)")
        end
        if !params.save_standardized_moments && !params.save_central_moments
            println("  Saving: Raw moments (M) only")
        end
    else
        println("\nSnapshots: DISABLED")
    end
    
    println("="^70)
end

# Export main functions
export parse_simulation_params, print_params_summary, get_default_params

