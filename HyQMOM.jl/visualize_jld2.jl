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
        @load filename snapshots grid params params_with_ic
        
        println("✓ Loaded successfully!")
        println("  Snapshots: $(length(snapshots))")
        println("  Time range: $(snapshots[1].t) to $(snapshots[end].t)")
        println("  Grid: $(params.Nx)×$(params.Ny)×$(params.Nz)")
        
        # Check what fields are available
        sample_snap = snapshots[1]
        available_fields = keys(sample_snap)
        println("  Available fields: $(collect(available_fields))")
        
        has_standardized = haskey(sample_snap, :S)
        if has_standardized
            println("  ✓ Standardized moments (S) available")
        else
            println("  ⚠ No standardized moments - run with --save-standardized-moments true for moment space view")
        end
        
        println("="^70)
        println("LAUNCHING INTERACTIVE 3D TIME-SERIES VIEWER")
        println("="^70)
        println("Controls:")
        println("  • Time slider: Navigate through snapshots")
        println("  • Play/Pause: Animate evolution")
        println("  • Quantity buttons: Switch between ρ, U, V, W, P")
        println("  • Iso Level slider: Adjust contour levels")
        println("  • Mouse: Rotate (drag), Zoom (scroll)")
        if has_standardized
            println("  • Moment space: Shows S110, S101, S011 correlations")
            println("  • Min |S| slider: Filter moment space points")
        end
        println("="^70)
        
        # Launch the viewer
        interactive_3d_timeseries(snapshots, grid, params_with_ic)
        
    catch e
        if isa(e, KeyError)
            println("ERROR: Invalid .jld2 file format")
            println("Expected fields: snapshots, grid, params, params_with_ic")
            println("This file may not be from HyQMOM examples/run_3d_custom_jets.jl")
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

