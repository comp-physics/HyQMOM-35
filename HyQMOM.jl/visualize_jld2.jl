#!/usr/bin/env julia

"""
Quick visualization script for HyQMOM .jld2 snapshot files

Usage:
    julia visualize_jld2.jl [filename.jld2]
    
If no filename provided, will interactively select from .jld2 files in current directory.

Examples:
    julia visualize_jld2.jl snapshots_crossing_Nx20_Ny20_Nz20_Kn1.00_Ma1.00_t0.0100.jld2
    julia visualize_jld2.jl
    
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

function find_jld2_files()
    files = filter(f -> endswith(f, ".jld2"), readdir("."))
    return files
end

function main()
    # Get filename from command line or find automatically
    if length(ARGS) >= 1
        filename = ARGS[1]
    else
        jld2_files = find_jld2_files()
        if isempty(jld2_files)
            println("ERROR: No .jld2 files found in current directory")
            println("Usage: julia visualize_jld2.jl [filename.jld2]")
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
            println("="^70)
        end
        
        # Launch the viewer with streaming file (outside the jldopen block)
        interactive_3d_timeseries_streaming(filename, grid, params_with_ic)
        
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

