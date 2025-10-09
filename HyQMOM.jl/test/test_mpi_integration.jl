"""
Integration Test - MPI Consistency

Tests that simulation results are identical between 1 and 2 MPI ranks.
This validates domain decomposition, halo exchange, and parallel correctness.

This test must be run in two stages:
1. First with 1 rank to generate reference data
2. Then with 2+ ranks to compare against reference

Usage:
    # Generate reference (1 rank)
    julia --project test/test_mpi_integration.jl
    
    # Compare (2 ranks)
    mpiexec -n 2 julia --project test/test_mpi_integration.jl
"""

using Test
using HyQMOM
using MPI

# Configuration
const Np_TEST = 20
const TMAX_TEST = 0.05
const MPI_TOL_ABS = 1e-10
const MPI_TOL_REL = 1e-8

@testset "MPI Integration Test" begin
    # Initialize MPI
    if !MPI.Initialized()
        MPI.Init()
    end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    if rank == 0
        println("\n=== Running MPI Integration Test ($(nprocs) rank(s)) ===")
    end
    
    # Run simulation
    results = run_simulation(
        Np = Np_TEST,
        tmax = TMAX_TEST,
        num_workers = nprocs,
        verbose = false,
        Kn = 1.0,
        Ma = 0.0,
        flag2D = 0,
        CFL = 0.5
    )
    
    if rank == 0
        M_final = results[:M]
        final_time = results[:final_time]
        time_steps = results[:time_steps]
        
        # Check for numerical issues
        @test !any(isnan, M_final)
        @test !any(isinf, M_final)
        
        # Basic sanity checks
        @test size(M_final) == (Np_TEST, Np_TEST, 35)
        @test final_time > 0
        @test time_steps > 0
        
        # Save reference data for 1 rank
        if nprocs == 1
            ref_file = joinpath(@__DIR__, "mpi_reference_$(Np_TEST).bin")
            println("Saving 1-rank reference to $(ref_file)")
            
            open(ref_file, "w") do io
                write(io, Int64(Np_TEST))
                write(io, Float64(final_time))
                write(io, Int64(time_steps))
                write(io, M_final)
            end
            
            println("Reference saved. Now run with multiple ranks to compare:")
            println("  mpiexec -n 2 julia --project test/test_mpi_integration.jl")
        else
            # Load and compare with 1-rank reference
            ref_file = joinpath(@__DIR__, "mpi_reference_$(Np_TEST).bin")
            
            if !isfile(ref_file)
                @warn "Reference file not found: $(ref_file)"
                @warn "First run with 1 rank to generate reference data:"
                @warn "  julia --project test/test_mpi_integration.jl"
                @test_skip "MPI consistency test (no reference data)"
            else
                println("Loading 1-rank reference from $(ref_file)")
                
                M_ref = nothing
                t_ref = 0.0
                steps_ref = 0
                
                open(ref_file, "r") do io
                    Np_saved = read(io, Int64)
                    t_ref = read(io, Float64)
                    steps_ref = read(io, Int64)
                    
                    M_ref = Array{Float64}(undef, Np_saved, Np_saved, 35)
                    read!(io, M_ref)
                end
                
                println("Comparing $(nprocs)-rank vs 1-rank results:")
                
                # Compare dimensions
                @test size(M_ref) == size(M_final)
                
                # Compare results
                diff = M_final .- M_ref
                abs_diff = abs.(diff)
                rel_diff = abs_diff ./ (abs.(M_ref) .+ 1e-30)
                
                max_abs_diff = maximum(abs_diff)
                max_rel_diff = maximum(rel_diff)
                mean_abs_diff = sum(abs_diff) / length(abs_diff)
                
                println("  Max absolute difference: $(max_abs_diff)")
                println("  Max relative difference: $(max_rel_diff)")
                println("  Mean absolute difference: $(mean_abs_diff)")
                
                # Find worst moment
                for k in 1:35
                    moment_max_abs = maximum(abs_diff[:, :, k])
                    if moment_max_abs > MPI_TOL_ABS * 10  # Report if > 10x tolerance
                        println("  WARNING: Moment $(k) has large diff: $(moment_max_abs)")
                    end
                end
                
                # Test tolerances
                @test max_abs_diff < MPI_TOL_ABS
                @test max_rel_diff < MPI_TOL_REL
                
                if max_abs_diff < MPI_TOL_ABS && max_rel_diff < MPI_TOL_REL
                    println("\n  OK MPI consistency test PASSED!")
                    println("     $(nprocs)-rank and 1-rank results match")
                else
                    println("\n  FAIL MPI consistency test FAILED!")
                    println("     Max abs diff: $(max_abs_diff) (tol: $(MPI_TOL_ABS))")
                    println("     Max rel diff: $(max_rel_diff) (tol: $(MPI_TOL_REL))")
                    
                    # Find location of max difference
                    max_idx = argmax(abs_diff)
                    max_loc = Tuple(max_idx)
                    println("     Location: $(max_loc)")
                    println("     1-rank:  $(M_ref[max_idx])")
                    println("     $(nprocs)-rank: $(M_final[max_idx])")
                end
            end
        end
    end
end

