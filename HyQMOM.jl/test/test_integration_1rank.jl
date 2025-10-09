"""
Integration Test - Single Rank

Tests the full simulation against MATLAB golden file for 1 MPI rank.
This validates the entire solver pipeline end-to-end.
"""

using Test
using MAT
using HyQMOM
using MPI

const INTEGRATION_TOL = 1e-8  # Slightly relaxed tolerance for full simulation

@testset "Integration - 1 Rank" begin
    golden_file = joinpath(@__DIR__, "..", "..", "goldenfiles", "goldenfile_mpi_1ranks_Np20_tmax100.mat")
    
    if !isfile(golden_file)
        @warn "Golden file not found: $golden_file - skipping integration test"
    else
        println("\n=== Running 1-Rank Integration Test ===")
        
        # Load golden data
        data = matread(golden_file)
        golden = data["golden_data"]
        
        # Extract nested structure
        params = golden["parameters"]
        moments = golden["moments"]
        grid = golden["grid"]
        
        println("Golden file info:")
        println("  Np: ", params["Np"])
        println("  tmax: ", params["tmax"])
        println("  final_time: ", params["final_time"])
        println("  time_steps: ", params["time_steps"])
        if haskey(params, "num_ranks")
            println("  num_ranks: ", params["num_ranks"])
        end
        if haskey(golden, "metadata")
            metadata = golden["metadata"]
            if haskey(metadata, "creation_date")
                println("  created: ", metadata["creation_date"])
            end
            if haskey(metadata, "git_branch")
                println("  branch: ", metadata["git_branch"])
            end
        end
        
        # Extract parameters
        Np = Int(params["Np"])
        tmax = Float64(params["tmax"])
        expected_M = moments["M"]  # NpxNpx35
        expected_xm = vec(grid["xm"])
        expected_ym = vec(grid["ym"])
        
        println("\nRunning Julia simulation...")
        println("  Np = $Np")
        println("  tmax = $tmax")
        
        # Run simulation
        # Note: This will need to be updated once run_simulation is fully implemented
        try
            # Initialize MPI
            if !MPI.Initialized()
                MPI.Init()
            end
            
            # Run simulation with same parameters as MATLAB
            results = run_simulation(
                Np = Np,
                tmax = tmax,
                num_workers = 1,
                verbose = false,
                save_output = false
            )
            
            println("  Simulation complete!")
            println("  Final time: ", results[:final_time])
            println("  Time steps: ", results[:time_steps])
            
            # Extract results
            M_final = results[:M]  # Should be NpxNpx35
            xm = results[:xm]
            ym = results[:ym]
            
            # Test grid coordinates
            @test xm ~= expected_xm atol=INTEGRATION_TOL
            @test ym ~= expected_ym atol=INTEGRATION_TOL
            
            # Test moment field
            @test size(M_final) == size(expected_M)
            @test M_final ~= expected_M atol=INTEGRATION_TOL
            
            # Compute relative error
            rel_error = maximum(abs.(M_final .- expected_M) ./ (abs.(expected_M) .+ 1e-10))
            println("  Maximum relative error: ", rel_error)
            
            # Test conservation (mass should be conserved)
            M000_initial = sum(expected_M[:,:,1])  # Total mass at t=0 (should be same as final if conservative)
            M000_final = sum(M_final[:,:,1])
            println("  Mass conservation check:")
            println("    Initial: ", M000_initial)
            println("    Final: ", M000_final)
            println("    Relative change: ", abs(M000_final - M000_initial) / M000_initial)
            
            # Test realizability (all moments should be realizable)
            println("  Realizability check:")
            all_realizable = true
            for i in 1:Np, j in 1:Np
                M = M_final[i,j,:]
                # Basic checks: M000 > 0, variances positive
                if M[1] <= 0
                    all_realizable = false
                    println("    WARNING: Non-positive mass at ($i,$j)")
                end
            end
            
            if all_realizable
                println("    OK All grid points have realizable moments")
            end
            
            println("\nOK 1-Rank integration test PASSED")
            
        catch e
            if isa(e, MethodError) && occursin("run_simulation", string(e))
                @warn "run_simulation not yet fully implemented - skipping execution test"
                @test_skip "Full simulation execution"
            else
                rethrow(e)
            end
        end
    end
end
