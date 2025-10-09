"""
MPI Golden File Consistency Tests

This test loads golden files generated with 1 rank and compares them against
simulations run with 2, 4, and 8 ranks to ensure MPI consistency.

Usage:
    mpiexec -n N julia --project=. test/test_mpi_goldenfiles.jl CONFIG

Where N is the number of ranks (2, 4, or 8) and CONFIG is "small" or "medium"
"""

using HyQMOM
using MPI
using Test
using Printf

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

# Get configuration from command line
config_name = length(ARGS) > 0 ? ARGS[1] : "small"

if rank == 0
    println("="^70)
    println("MPI GOLDEN FILE TEST: $(nprocs) ranks, config=$(config_name)")
    println("="^70)
    println()
end

# Load golden file
golden_dir = joinpath(@__DIR__, "goldenfiles")
golden_file = joinpath(golden_dir, "mpi_1rank_$(config_name).bin")

if !isfile(golden_file)
    if rank == 0
        println("ERROR: Golden file not found: $(golden_file)")
        println("Please run: mpiexec -n 1 julia --project=. test/create_mpi_goldenfiles.jl")
    end
    MPI.Finalize()
    exit(1)
end

# Load golden data
golden_nranks, golden_Np, golden_tmax, golden_tfinal, golden_steps, M_golden = 
    open(golden_file, "r") do io
        nranks = read(io, Int64)
        Np = read(io, Int64)
        tmax = read(io, Float64)
        tfinal = read(io, Float64)
        steps = read(io, Int64)
        M = Array{Float64}(undef, Np, Np, 35)
        read!(io, M)
        (nranks, Np, tmax, tfinal, steps, M)
    end

if rank == 0
    println("Loaded golden file:")
    println("  Ranks: $(golden_nranks)")
    println("  Grid: $(golden_Np)x$(golden_Np)")
    println("  tmax: $(golden_tmax)")
    println("  Final time: $(golden_tfinal)")
    println("  Steps: $(golden_steps)")
    println()
end

# Run simulation with current number of ranks
if rank == 0
    println("Running $(nprocs)-rank simulation...")
end

results = run_simulation(
    Np = golden_Np,
    tmax = golden_tmax,
    num_workers = nprocs,
    verbose = false,
    Kn = 1.0,
    Ma = 0.0,
    flag2D = 0,
    CFL = 0.5
)

if rank == 0
    M_test = results[:M]
    t_final = results[:final_time]
    steps = results[:time_steps]
    
    println("Simulation complete:")
    println("  Final time: $(t_final)")
    println("  Steps: $(steps)")
    println()
    
    # Compare results
    @testset "MPI Golden File: $(nprocs) ranks vs 1 rank ($(config_name))" begin
        # Check dimensions
        @test size(M_test) == size(M_golden)
        
        # Check final time (should be very close)
        @test t_final ≈ golden_tfinal atol=1e-10
        
        # Check moments
        diff_abs = abs.(M_test .- M_golden)
        max_abs_diff = maximum(diff_abs)
        mean_abs_diff = sum(diff_abs) / length(diff_abs)
        
        # For relative difference, avoid division by small numbers
        mask = abs.(M_golden) .> 1e-10
        rel_diff = zeros(size(M_golden))
        rel_diff[mask] = abs.((M_test[mask] .- M_golden[mask]) ./ M_golden[mask])
        max_rel_diff = maximum(rel_diff)
        
        println("Comparison Results:")
        println("  Max absolute difference: $(max_abs_diff)")
        println("  Mean absolute difference: $(mean_abs_diff)")
        println("  Max relative difference: $(max_rel_diff)")
        println()
        
        # Test with tight tolerance (machine epsilon)
        @test max_abs_diff < 1e-13
        @test max_rel_diff < 1e-11
        
        if max_abs_diff < 1e-13 && max_rel_diff < 1e-11
            println("  ✓ MPI consistency test PASSED!")
            println("    $(nprocs)-rank and 1-rank results match within machine epsilon")
        else
            println("  ✗ MPI consistency test FAILED!")
            println("    Differences exceed expected tolerance")
        end
    end
end

MPI.Finalize()

