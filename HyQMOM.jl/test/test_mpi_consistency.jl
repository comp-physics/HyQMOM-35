"""
MPI Consistency Test

Tests that simulation results are consistent between 1 and 2 MPI ranks.
This validates that the domain decomposition and halo exchange work correctly.

This test should be run with:
    julia --project test/test_mpi_consistency.jl         # 1 rank
    mpiexec -n 2 julia --project test/test_mpi_consistency.jl  # 2 ranks
"""

using Test
using HyQMOM
using MPI
using Printf

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

# Only rank 0 prints
function print_rank0(msg)
    if rank == 0
        println(msg)
    end
end

# Test configuration
const TOLERANCE_ABS = 1e-10  # Very tight tolerance for deterministic results
const TOLERANCE_REL = 1e-8

# Simple parameters for quick test
Np = 20
tmax = 0.05  # Short simulation for fast testing
Kn = 1.0
Ma = 0.0
flag2D = 0
CFL = 0.5

print_rank0("="^70)
print_rank0("MPI CONSISTENCY TEST: $(nprocs) ranks")
print_rank0("="^70)
print_rank0("")
print_rank0("This script should be run twice:")
print_rank0("  1. julia --project test/test_mpi_consistency.jl")
print_rank0("  2. mpiexec -n 2 julia --project test/test_mpi_consistency.jl")
print_rank0("")
print_rank0("Configuration:")
print_rank0("  Np = $Np")
print_rank0("  tmax = $tmax")
print_rank0("  Kn = $Kn, Ma = $Ma")
print_rank0("  CFL = $CFL")
print_rank0("="^70)

# Run simulation
print_rank0("\n[RUN] Running simulation with $nprocs rank(s)...")

results = run_simulation(
    Np = Np,
    tmax = tmax,
    num_workers = nprocs,
    verbose = false,
    save_output = false,
    Kn = Kn,
    Ma = Ma,
    flag2D = flag2D,
    CFL = CFL
)

if rank == 0
    M_final = results[:M]
    final_time = results[:final_time]
    time_steps = results[:time_steps]
    
    println("  OK Simulation complete!")
    println("     Final time: $(final_time)")
    println("     Time steps: $(time_steps)")
    println("     Result size: ", size(M_final))
    
    # Save results to file
    output_file = "mpi_results_$(nprocs)ranks_Np$(Np).bin"
    println("\n[SAVE] Saving results to $(output_file)...")
    
    open(output_file, "w") do io
        # Write metadata
        write(io, Int64(nprocs))
        write(io, Int64(Np))
        write(io, Float64(final_time))
        write(io, Int64(time_steps))
        
        # Write moment array
        write(io, M_final)
    end
    
    println("  OK Saved successfully")
    
    # If we have results from both 1 and 2 ranks, compare them
    file_1rank = "mpi_results_1ranks_Np$(Np).bin"
    file_2ranks = "mpi_results_2ranks_Np$(Np).bin"
    
    if nprocs == 2 && isfile(file_1rank)
        println("\n" * "="^70)
        println("COMPARISON: 1 rank vs 2 ranks")
        println("="^70)
        
        # Load 1-rank results
        println("\n[LOAD] Loading 1-rank results from $(file_1rank)...")
        M_1rank = nothing
        t_1rank = 0.0
        steps_1rank = 0
        
        open(file_1rank, "r") do io
            nprocs_saved = read(io, Int64)
            Np_saved = read(io, Int64)
            t_1rank = read(io, Float64)
            steps_1rank = read(io, Int64)
            
            M_1rank = Array{Float64}(undef, Np_saved, Np_saved, 35)
            read!(io, M_1rank)
        end
        
        println("  OK Loaded 1-rank: t=$(t_1rank), steps=$(steps_1rank)")
        
        # Compare dimensions
        println("\n[CHECK] Comparing results...")
        if size(M_1rank) != size(M_final)
            println("  FAIL DIMENSION MISMATCH!")
            println("     1 rank:  ", size(M_1rank))
            println("     2 ranks: ", size(M_final))
            @test false  # Fail the test
        else
            println("  OK Dimensions match: ", size(M_final))
        end
        
        # Compare time and steps
        if abs(t_1rank - final_time) > 1e-12
            println("  WARNING  Final time differs slightly:")
            println("     1 rank:  $(t_1rank)")
            println("     2 ranks: $(final_time)")
            println("     diff:    $(abs(t_1rank - final_time))")
        else
            println("  OK Final time matches: t=$(final_time)")
        end
        
        if steps_1rank != time_steps
            println("  WARNING  Time step count differs:")
            println("     1 rank:  $(steps_1rank)")
            println("     2 ranks: $(time_steps)")
        else
            println("  OK Time step count matches: $(time_steps) steps")
        end
        
        # Compute differences
        diff = M_final .- M_1rank
        abs_diff = abs.(diff)
        rel_diff = abs_diff ./ (abs.(M_1rank) .+ 1e-30)
        
        max_abs_diff = maximum(abs_diff)
        max_rel_diff = maximum(rel_diff)
        mean_abs_diff = sum(abs_diff) / length(abs_diff)
        mean_rel_diff = sum(rel_diff) / length(rel_diff)
        
        # Find location of maximum difference
        max_diff_idx = argmax(abs_diff)
        max_diff_loc = Tuple(max_diff_idx)
        
        println("\n  Absolute Differences:")
        println("    Max:  ", @sprintf("%.6e", max_abs_diff), " at ", max_diff_loc)
        println("    Mean: ", @sprintf("%.6e", mean_abs_diff))
        
        println("\n  Relative Differences:")
        println("    Max:  ", @sprintf("%.6e", max_rel_diff))
        println("    Mean: ", @sprintf("%.6e", mean_rel_diff))
        
        # Check for NaN or Inf
        nan_count_1 = sum(isnan, M_1rank)
        inf_count_1 = sum(isinf, M_1rank)
        nan_count_2 = sum(isnan, M_final)
        inf_count_2 = sum(isinf, M_final)
        
        if nan_count_1 > 0 || inf_count_1 > 0 || nan_count_2 > 0 || inf_count_2 > 0
            println("\n  FAIL NUMERICAL ISSUES:")
            println("     1 rank  - NaN: $(nan_count_1), Inf: $(inf_count_1)")
            println("     2 ranks - NaN: $(nan_count_2), Inf: $(inf_count_2)")
            @test false
        else
            println("\n  OK No NaN or Inf values in either result")
        end
        
        # Detailed statistics by moment
        println("\n  Detailed Differences by Moment Component:")
        println("  " * "-"^66)
        nx, ny, nmom = size(M_final)
        
        max_moment_diff = 0.0
        worst_moment = 0
        
        for k in 1:nmom
            moment_abs_diff = abs_diff[:, :, k]
            moment_rel_diff = rel_diff[:, :, k]
            
            max_abs = maximum(moment_abs_diff)
            max_rel = maximum(moment_rel_diff)
            mean_abs = sum(moment_abs_diff) / length(moment_abs_diff)
            
            if max_abs > max_moment_diff
                max_moment_diff = max_abs
                worst_moment = k
            end
            
            # Print all moments
            @printf("    Moment %2d: max_abs=%.6e, max_rel=%.6e, mean_abs=%.6e\n",
                    k, max_abs, max_rel, mean_abs)
        end
        
        # Overall pass/fail
        println("\n" * "="^70)
        if max_abs_diff < TOLERANCE_ABS && max_rel_diff < TOLERANCE_REL
            println("OK MPI CONSISTENCY TEST PASSED!")
            println("   1-rank and 2-rank results match within tolerance")
            println("   Absolute tolerance: ", @sprintf("%.6e", TOLERANCE_ABS))
            println("   Relative tolerance: ", @sprintf("%.6e", TOLERANCE_REL))
            @test true
        else
            println("FAIL MPI CONSISTENCY TEST FAILED!")
            println("   Differences exceed tolerance:")
            println("   Max abs diff: ", @sprintf("%.6e", max_abs_diff), 
                    " (tolerance: ", @sprintf("%.6e", TOLERANCE_ABS), ")")
            println("   Max rel diff: ", @sprintf("%.6e", max_rel_diff),
                    " (tolerance: ", @sprintf("%.6e", TOLERANCE_REL), ")")
            println("   Worst moment: #$(worst_moment)")
            
            # Show details at max diff location
            println("\n   At max diff location ", max_diff_loc, ":")
            println("     1 rank:  ", @sprintf("%.10e", M_1rank[max_diff_idx]))
            println("     2 ranks: ", @sprintf("%.10e", M_final[max_diff_idx]))
            println("     Diff:    ", @sprintf("%.10e", diff[max_diff_idx]))
            
            @test false
        end
        println("="^70)
        
    elseif nprocs == 1
        println("\n" * "="^70)
        println("NEXT STEP: Run with 2 ranks to compare:")
        println("  mpiexec -n 2 julia --project test/test_mpi_consistency.jl")
        println("="^70)
    end
end

# Finalize MPI
MPI.Finalize()

