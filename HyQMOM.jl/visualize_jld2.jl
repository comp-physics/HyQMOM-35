#!/usr/bin/env julia

"""
Quick visualization script for HyQMOM .jld2 snapshot files

Usage:
    julia visualize_jld2.jl [filename.jld2] [--snapshot first|last|N]
    
If no filename provided, will interactively select from .jld2 files in current directory.

Options:
    --snapshot first    Show only the first snapshot
    --snapshot last     Show only the last snapshot
    --snapshot N        Show only snapshot number N (1-indexed)
    (default: show all snapshots with time slider)

Examples:
    julia visualize_jld2.jl snapshots.jld2                    # All snapshots
    julia visualize_jld2.jl snapshots.jld2 --snapshot first   # First only
    julia visualize_jld2.jl snapshots.jld2 --snapshot last    # Last only
    julia visualize_jld2.jl snapshots.jld2 --snapshot 5       # Snapshot #5 only
    julia visualize_jld2.jl                                   # Interactive file selection
    
Generate snapshot files using:
    julia examples/run_3d_crossing_jets.jl --save-standardized-moments true
    julia examples/run_3d_jets_timeseries.jl --save-standardized-moments true
"""

using Pkg
Pkg.activate(@__DIR__)

# Check if GLMakie is installed
try
    using GLMakie
catch e
    println("ERROR: GLMakie not installed. Installing now...")
    Pkg.add("GLMakie")
    using GLMakie
end

using HyQMOM, JLD2
using Printf

function show_help()
    println("""
    HyQMOM JLD2 Visualization Tool
    ==============================
    
    USAGE:
        julia visualize_jld2.jl [OPTIONS] [FILENAME]
    
    ARGUMENTS:
        FILENAME                JLD2 snapshot file to visualize
                               If omitted, will prompt for file selection
    
    OPTIONS:
        --help, -h             Show this help message and exit
        
        --snapshot MODE        Control which snapshot(s) to display:
          first                Show only the first snapshot
          last                 Show only the last snapshot  
          N                    Show only snapshot number N (1-indexed)
                               (default: show all snapshots with time slider)
    
    VISUALIZATION FEATURES:
        * Interactive 3D isosurface rendering
        * Time-series animation with play/pause controls
        * Multiple quantities: Density, Velocity (U/V/W), Pressure
        * Moment space visualization (if standardized moments available)
        * Mouse controls: Rotate (drag), Zoom (scroll)
        * Adjustable iso-levels and transparency
    
    EXAMPLES:
        # Interactive file selection with all snapshots
        julia visualize_jld2.jl
        
        # Visualize all snapshots from specific file
        julia visualize_jld2.jl snapshots.jld2
        
        # Show only the first snapshot
        julia visualize_jld2.jl snapshots.jld2 --snapshot first
        
        # Show only the last snapshot  
        julia visualize_jld2.jl snapshots.jld2 --snapshot last
        
        # Show specific snapshot (e.g., snapshot #5)
        julia visualize_jld2.jl snapshots.jld2 --snapshot 5
    
    GENERATING SNAPSHOT FILES:
        Run examples with --save-standardized-moments flag:
        julia examples/run_3d_crossing_jets.jl --save-standardized-moments true
        julia examples/run_3d_jets_timeseries.jl --save-standardized-moments true
    
    REQUIREMENTS:
        * GLMakie package (installed automatically if missing)
        * JLD2 snapshot file from HyQMOM simulation
    """)
end

function find_jld2_files()
    files = filter(f -> endswith(f, ".jld2"), readdir("."))
    return files
end

function parse_args()
    filename = nothing
    snapshot_mode = :all  # :all, :first, :last, or specific number
    snapshot_number = nothing
    
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--help" || arg == "-h"
            show_help()
            exit(0)
        elseif arg == "--snapshot"
            if i + 1 > length(ARGS)
                error("--snapshot requires an argument (first|last|N)")
            end
            snapshot_arg = ARGS[i + 1]
            if snapshot_arg == "first"
                snapshot_mode = :first
            elseif snapshot_arg == "last"
                snapshot_mode = :last
            else
                # Try to parse as number
                try
                    snapshot_number = parse(Int, snapshot_arg)
                    snapshot_mode = :specific
                catch
                    error("Invalid --snapshot argument: $snapshot_arg (expected first|last|N)")
                end
            end
            i += 2
        elseif startswith(arg, "--")
            println("ERROR: Unknown option: $arg")
            println("Run 'julia visualize_jld2.jl --help' for usage information")
            exit(1)
        else
            # Assume it's the filename
            filename = arg
            i += 1
        end
    end
    
    return filename, snapshot_mode, snapshot_number
end

function main()
    # Parse command line arguments
    filename, snapshot_mode, snapshot_number = parse_args()
    
    # Get filename from command line or find automatically
    if filename === nothing
        jld2_files = find_jld2_files()
        if isempty(jld2_files)
            println("ERROR: No .jld2 files found in current directory")
            println("Usage: julia visualize_jld2.jl [filename.jld2] [--snapshot first|last|N]")
            return
        elseif length(jld2_files) == 1
            filename = jld2_files[1]
            println("Found single .jld2 file: $filename")
        else
            println("Multiple .jld2 files found:")
            for (i, f) in enumerate(jld2_files)
                println("  $i: $f")
            end
            print("Select file number (1-$(length(jld2_files))): ")
            choice = parse(Int, readline())
            filename = jld2_files[choice]
        end
    end
    
    # Check if file exists
    if !isfile(filename)
        println("ERROR: File not found: $filename")
        return
    end
    
    println("="^70)
    println("LOADING: $filename")
    if snapshot_mode != :all
        if snapshot_mode == :first
            println("MODE: Showing first snapshot only")
        elseif snapshot_mode == :last
            println("MODE: Showing last snapshot only")
        else
            println("MODE: Showing snapshot #$snapshot_number only")
        end
    end
    println("="^70)
    
    # Load the data
    try
        # Variables to store outside jldopen block
        local grid, params_with_ic, has_standardized
        
        jldopen(filename, "r") do f
            # Check format and extract metadata
            if haskey(f, "meta/n_snapshots")
                # Streaming format (current)
                n_snapshots = f["meta/n_snapshots"]
                grid = f["grid"]
                params = f["meta/params"]
                params_with_ic = params
                
                println("[OK] Loaded successfully!")
                println("  Snapshots: $n_snapshots")
                
                # Load first and last snapshot to get time range
                snap_keys = sort!(collect(keys(f["snapshots"])))
                first_snap = f["snapshots/$(snap_keys[1])"]
                last_snap = f["snapshots/$(snap_keys[end])"]
                
                println("  Time range: $(first_snap["t"]) to $(last_snap["t"])")
                println("  Grid: $(params.Nx)x$(params.Ny)x$(params.Nz)")
                
                # Validate snapshot selection
                if snapshot_mode == :specific
                    if snapshot_number < 1 || snapshot_number > n_snapshots
                        println("ERROR: Snapshot number $snapshot_number out of range (1-$n_snapshots)")
                        return
                    end
                end
                
                # Check what fields are available
                has_standardized = haskey(first_snap, "S")
                if has_standardized
                    println("  [OK] Standardized moments (S) available")
                else
                    println("  [WARNING] No standardized moments - run with --save-standardized-moments true for moment space view")
                end
            else
                # Legacy format - convert to streaming by saving
                println("[WARNING] Legacy format detected - converting to streaming format...")
                @load filename snapshots grid params params_with_ic
                
                # Save in streaming format
                new_filename = replace(filename, r"\.jld2$" => "_streaming.jld2")
                jldopen(new_filename, "w") do f_new
                    f_new["meta/params"] = params
                    f_new["meta/snapshot_interval"] = 1  # Unknown, use 1
                    f_new["meta/n_snapshots"] = length(snapshots)
                    f_new["grid"] = grid
                    
                    for (i, snap) in enumerate(snapshots)
                        snap_key = lpad(i, 6, '0')
                        f_new["snapshots/$snap_key/M"] = snap.M
                        f_new["snapshots/$snap_key/t"] = snap.t
                        f_new["snapshots/$snap_key/step"] = snap.step
                        if haskey(snap, :S)
                            f_new["snapshots/$snap_key/S"] = snap.S
                        end
                        if haskey(snap, :C)
                            f_new["snapshots/$snap_key/C"] = snap.C
                        end
                    end
                end
                
                println("[OK] Converted to streaming format: $new_filename")
                filename = new_filename
                
                params_with_ic = params
                has_standardized = haskey(snapshots[1], :S)
            end
            
            println("="^70)
            if snapshot_mode == :all
                println("LAUNCHING INTERACTIVE 3D TIME-SERIES VIEWER")
                println("="^70)
                println("Controls:")
                println("  * Time slider: Navigate through snapshots (loaded on-demand)")
                println("  * Play/Pause: Animate evolution")
                println("  * Quantity buttons: Switch between rho, U, V, W, P")
                println("  * Iso Level slider: Adjust contour levels")
                println("  * Mouse: Rotate (drag), Zoom (scroll)")
                if has_standardized
                    println("  * Moment space: Shows S110, S101, S011 correlations")
                    println("  * Min |S| slider: Filter moment space points")
                end
            else
                println("LAUNCHING SINGLE SNAPSHOT VIEWER")
                println("="^70)
                println("Controls:")
                println("  * Quantity buttons: Switch between rho, U, V, W, P")
                println("  * Iso Level slider: Adjust contour levels")
                println("  * Mouse: Rotate (drag), Zoom (scroll)")
                if has_standardized
                    println("  * Moment space: Shows S110, S101, S011 correlations")
                    println("  * Min |S| slider: Filter moment space points")
                end
            end
            println("="^70)
        end
        
        # Launch the appropriate viewer (outside the jldopen block)
        if snapshot_mode == :all
            # Show all snapshots with time slider
            interactive_3d_timeseries_streaming(filename, grid, params_with_ic)
        else
            # Show single snapshot
            interactive_3d_timeseries_streaming(filename, grid, params_with_ic, 
                                               snapshot_mode=snapshot_mode,
                                               snapshot_number=snapshot_number)
        end
        
    catch e
        if isa(e, KeyError)
            println("ERROR: Invalid .jld2 file format")
            println("Expected fields: meta/n_snapshots, grid, meta/params OR snapshots, grid, params")
            println("This file may not be from HyQMOM examples (e.g., run_3d_crossing_jets.jl, run_3d_jets_timeseries.jl)")
        else
            println("ERROR loading file: $e")
            rethrow()
        end
        return
    end
end

# Run if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

