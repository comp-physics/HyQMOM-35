#!/usr/bin/env julia
"""
Test Julia implementation against MATLAB golden files with exact parameter matching.

This test ensures that the Julia code produces identical results to MATLAB by:
1. Using the exact same parameters as MATLAB (Kn=1.0, Ma=0.0, CFL=0.5, etc.)
2. Running with the same grid size and time
3. Comparing against the MATLAB golden files

Run with:
    julia --project test/test_matlab_golden.jl
    mpiexec -n 1 julia --project test/test_matlab_golden.jl
"""

using Test
using Printf
using MAT
using MPI

# Add parent directory to load path
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using RodneyHQMOM

function load_matlab_golden(golden_file)
    """Load MATLAB golden file and extract data."""
    if !isfile(golden_file)
        error("Golden file not found: $golden_file")
    end
    
    println("📥 Loading MATLAB golden file: $golden_file")
    matlab_data = matread(golden_file)
    golden_data = matlab_data["golden_data"]
    
    # Extract from nested structure
    M_matlab = golden_data["moments"]["M"]
    Np = Int(golden_data["parameters"]["Np"])
    tmax = golden_data["parameters"]["tmax"]
    final_time = golden_data["parameters"]["final_time"]
    time_steps = Int(golden_data["parameters"]["time_steps"])
    num_ranks = Int(golden_data["parameters"]["num_workers"])
    
    println("  ✓ Loaded: Np=$Np, tmax=$tmax, ranks=$num_ranks")
    println("  ✓ MATLAB final time: $final_time ($(time_steps) steps)")
    println("  ✓ Moment array shape: ", size(M_matlab))
    
    return (M=M_matlab, Np=Np, tmax=tmax, final_time=final_time, 
            time_steps=time_steps, num_ranks=num_ranks)
end

function run_julia_simulation_exact_matlab_params(Np, tmax, num_workers=1)
    """
    Run Julia simulation with EXACT MATLAB parameters.
    
    MATLAB uses these parameters (from parse_main_args.m):
    - Kn = 1.0
    - Ma = 0.0
    - CFL = 0.5
    - rhol = 1.0
    - rhor = 0.01
    - T = 1.0
    - r110, r101, r011 = 0.0
    - flag2D = 0 (3D)
    """
    
    println("\n🚀 Running Julia simulation with EXACT MATLAB parameters:")
    println("  Np = $Np")
    println("  tmax = $tmax")
    println("  Kn = 1.0  (MATLAB default)")
    println("  Ma = 0.0  (MATLAB default)")
    println("  CFL = 0.5 (MATLAB default)")
    println("  rhol = 1.0, rhor = 0.01 (MATLAB IC)")
    println("  T = 1.0")
    println("  flag2D = 0 (3D)")
    
    # Setup parameters EXACTLY as MATLAB does
    dx = 1.0 / Np
    dy = 1.0 / Np
    Nmom = 35
    nnmax = 20000000  # MATLAB uses 2e7
    Kn = 1.0
    dtmax = Kn  # MATLAB sets dtmax = Kn
    
    # Initial condition parameters - EXACTLY as MATLAB
    rhol = 1.0
    rhor = 0.01  # MATLAB uses 0.01, not 0.5!
    T = 1.0
    r110 = 0.0
    r101 = 0.0
    r011 = 0.0
    
    # Physical parameters - EXACTLY as MATLAB
    Ma = 0.0   # MATLAB uses 0.0
    flag2D = 0
    CFL = 0.5  # MATLAB uses 0.5, not 0.9!
    
    # Diagnostic parameters
    symmetry_check_interval = 10
    
    # Package parameters into named tuple
    params = (
        Np=Np, 
        tmax=tmax, 
        Kn=Kn, 
        Ma=Ma, 
        flag2D=flag2D, 
        CFL=CFL,
        dx=dx, 
        dy=dy, 
        Nmom=Nmom, 
        nnmax=nnmax, 
        dtmax=dtmax,
        rhol=rhol, 
        rhor=rhor, 
        T=T, 
        r110=r110, 
        r101=r101, 
        r011=r011,
        symmetry_check_interval=symmetry_check_interval,
        enable_memory_tracking=false
    )
    
    # Run simulation
    start_time = time()
    M_final, final_time, time_steps, grid_out = RodneyHQMOM.simulation_runner(params)
    elapsed = time() - start_time
    
    println("  ✓ Julia simulation complete in $(round(elapsed, digits=2))s")
    println("  ✓ Julia final time: $final_time ($(time_steps) steps)")
    if M_final !== nothing
        println("  ✓ Moment array shape: ", size(M_final))
    end
    
    return (M=M_final, final_time=final_time, time_steps=time_steps, grid=grid_out)
end

function compare_results(matlab_result, julia_result, tolerance_abs=1e-8, tolerance_rel=1e-6)
    """Compare MATLAB and Julia results with detailed diagnostics."""
    
    println("\n" * "="^70)
    println("📊 COMPARISON RESULTS")
    println("="^70)
    
    M_matlab = matlab_result.M
    M_julia = julia_result.M
    
    if M_julia === nothing
        println("  ⚠️  Julia result is nothing (not on rank 0)")
        return false
    end
    
    # Check dimensions
    println("\n1. Dimension Check:")
    if size(M_matlab) != size(M_julia)
        println("  ❌ DIMENSION MISMATCH!")
        println("     MATLAB: ", size(M_matlab))
        println("     Julia:  ", size(M_julia))
        return false
    else
        println("  ✓ Dimensions match: ", size(M_julia))
    end
    
    # Check time
    println("\n2. Time Check:")
    time_diff = abs(matlab_result.final_time - julia_result.final_time)
    if time_diff > 1e-10
        println("  ⚠️  Time mismatch:")
        println("     MATLAB: $(matlab_result.final_time) ($(matlab_result.time_steps) steps)")
        println("     Julia:  $(julia_result.final_time) ($(julia_result.time_steps) steps)")
        println("     Diff:   $time_diff")
    else
        println("  ✓ Time matches: $(julia_result.final_time)")
        println("     MATLAB steps: $(matlab_result.time_steps)")
        println("     Julia steps:  $(julia_result.time_steps)")
    end
    
    # Check for NaN or Inf
    println("\n3. Numerical Health Check:")
    nan_count = sum(isnan, M_julia)
    inf_count = sum(isinf, M_julia)
    
    if nan_count > 0 || inf_count > 0
        println("  ❌ NUMERICAL ISSUES:")
        println("     NaN count: $nan_count")
        println("     Inf count: $inf_count")
        return false
    else
        println("  ✓ No NaN or Inf values")
    end
    
    # Compute differences
    println("\n4. Value Comparison:")
    diff = M_julia .- M_matlab
    abs_diff = abs.(diff)
    rel_diff = abs_diff ./ (abs.(M_matlab) .+ 1e-30)
    
    max_abs_diff = maximum(abs_diff)
    max_rel_diff = maximum(rel_diff)
    mean_abs_diff = sum(abs_diff) / length(abs_diff)
    mean_rel_diff = sum(rel_diff) / length(rel_diff)
    
    # Find location of maximum difference
    max_diff_idx = argmax(abs_diff)
    max_diff_loc = Tuple(max_diff_idx)
    
    println("\n  Absolute Differences:")
    @printf("    Max:  %.6e at %s\n", max_abs_diff, max_diff_loc)
    @printf("    Mean: %.6e\n", mean_abs_diff)
    
    println("\n  Relative Differences:")
    @printf("    Max:  %.6e\n", max_rel_diff)
    @printf("    Mean: %.6e\n", mean_rel_diff)
    
    println("\n  At max diff location $max_diff_loc:")
    @printf("    MATLAB: %.10e\n", M_matlab[max_diff_idx])
    @printf("    Julia:  %.10e\n", M_julia[max_diff_idx])
    @printf("    Diff:   %.10e\n", diff[max_diff_idx])
    
    # Detailed statistics by moment
    println("\n5. Moment-by-Moment Analysis:")
    println("-"^70)
    nx, ny, nmom = size(M_julia)
    
    moment_names = ["ρ", "ρu", "ρv", "ρw", "ρE", "M₁₁₀", "M₁₀₁", "M₀₁₁", 
                    "M₂₀₀", "M₀₂₀", "M₀₀₂"]
    
    for k in 1:min(11, nmom)
        moment_abs_diff = abs_diff[:, :, k]
        moment_rel_diff = rel_diff[:, :, k]
        
        max_abs = maximum(moment_abs_diff)
        max_rel = maximum(moment_rel_diff)
        mean_abs = sum(moment_abs_diff) / length(moment_abs_diff)
        
        moment_name = k <= length(moment_names) ? moment_names[k] : "M[$k]"
        
        @printf("  %-6s: max_abs=%.6e, max_rel=%.6e, mean_abs=%.6e", 
                moment_name, max_abs, max_rel, mean_abs)
        
        # Flag if this moment has issues
        if max_abs > tolerance_abs || max_rel > tolerance_rel
            print("  ⚠️")
        end
        println()
    end
    
    # Overall pass/fail
    println("\n" * "="^70)
    println("6. Final Verdict:")
    println("-"^70)
    @printf("  Tolerance (absolute): %.6e\n", tolerance_abs)
    @printf("  Tolerance (relative): %.6e\n", tolerance_rel)
    @printf("  Max absolute diff:    %.6e\n", max_abs_diff)
    @printf("  Max relative diff:    %.6e\n", max_rel_diff)
    println()
    
    passed = (max_abs_diff < tolerance_abs && max_rel_diff < tolerance_rel)
    
    if passed
        println("  ✅ TEST PASSED!")
        println("  Julia matches MATLAB within tolerance.")
    else
        println("  ❌ TEST FAILED!")
        println("  Differences exceed tolerance thresholds.")
        
        # Additional diagnostics for failure
        exceed_abs = sum(abs_diff .> tolerance_abs)
        exceed_rel = sum(rel_diff .> tolerance_rel)
        total_points = length(abs_diff)
        
        @printf("\n  Points exceeding absolute tolerance: %d / %d (%.2f%%)\n", 
                exceed_abs, total_points, 100*exceed_abs/total_points)
        @printf("  Points exceeding relative tolerance: %d / %d (%.2f%%)\n",
                exceed_rel, total_points, 100*exceed_rel/total_points)
    end
    
    println("="^70)
    
    return passed
end

function main()
    println("="^70)
    println("MATLAB GOLDEN FILE TEST - Julia vs MATLAB Exact Comparison")
    println("="^70)
    println()
    println("This test runs Julia with EXACT MATLAB parameters:")
    println("  Kn=1.0, Ma=0.0, CFL=0.5, rhol=1.0, rhor=0.01")
    println()
    
    # Initialize MPI
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    if rank == 0
        println("MPI Configuration:")
        println("  Ranks: $nprocs")
        println()
    end
    
    # Test configuration
    golden_file = "../goldenfiles/goldenfile_mpi_1ranks_Np20_tmax100.mat"
    
    # Load MATLAB golden file
    matlab_result = nothing
    if rank == 0
        try
            matlab_result = load_matlab_golden(golden_file)
        catch e
            println("❌ Error loading golden file: $e")
            MPI.Finalize()
            return 1
        end
    end
    
    # Broadcast parameters to all ranks
    Np = 20
    tmax = 0.1
    
    if rank == 0
        println("\n" * "="^70)
        println("Running Julia simulation with 1 MPI rank")
        println("="^70)
    end
    
    # Run Julia simulation with exact MATLAB parameters
    julia_result = run_julia_simulation_exact_matlab_params(Np, tmax, 1)
    
    # Compare results (only on rank 0)
    test_passed = false
    if rank == 0 && matlab_result !== nothing
        test_passed = compare_results(matlab_result, julia_result)
        
        println("\n" * "="^70)
        if test_passed
            println("🎉 OVERALL TEST RESULT: PASSED")
        else
            println("❌ OVERALL TEST RESULT: FAILED")
            println("\nPossible issues:")
            println("  1. Different random number generation")
            println("  2. Different order of operations in collision operator")
            println("  3. Different realizability enforcement")
            println("  4. Different flux computation")
            println("\nNext steps:")
            println("  - Check collision operator implementation")
            println("  - Verify realizability functions")
            println("  - Compare intermediate values during time stepping")
        end
        println("="^70)
    end
    
    MPI.Finalize()
    
    return test_passed ? 0 : 1
end

# Run if this is the main script
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end

