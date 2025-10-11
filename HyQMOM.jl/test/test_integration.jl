"""
Full Simulation Integration Tests

Tests the complete simulation pipeline against MATLAB golden files.
This consolidates what was previously in:
- test_matlab_golden.jl
- test_matlab_golden_simple.jl  
- test_integration_1rank.jl

Usage:
    # Within Pkg.test() framework
    julia --project=. -e 'using Pkg; Pkg.test()'
    
    # Standalone with MPI
    julia --project=. test/test_integration.jl
    mpiexec -n 1 julia --project=. test/test_integration.jl
    
    # Skip with environment variable
    TEST_INTEGRATION=false julia --project=. -e 'using Pkg; Pkg.test()'
"""

using Test
using Printf
using MAT

# Add HyQMOM to load path and import it
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using HyQMOM

# Check if MPI is available
const HAS_MPI = try
    using MPI
    true
catch
    false
end

# Configuration
const INTEGRATION_TOL_ABS = 1e-8
const INTEGRATION_TOL_REL = 1e-6
const GOLDEN_FILE = joinpath(@__DIR__, "..", "..", "goldenfiles", 
                              "goldenfile_mpi_1ranks_Np20_tmax100.mat")
const RUN_INTEGRATION = get(ENV, "TEST_INTEGRATION", "true") != "false"

# Check if we're running standalone or within Test framework
const STANDALONE = abspath(PROGRAM_FILE) == @__FILE__

# Helper function to initialize MPI if needed
function ensure_mpi_initialized()
    if !HAS_MPI
        return (rank=0, nprocs=1, comm=nothing)
    end
    
    if !MPI.Initialized()
        MPI.Init()
    end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    return (rank=rank, nprocs=nprocs, comm=comm)
end

function load_matlab_golden(golden_file)
    """Load MATLAB golden file and extract data."""
    if !isfile(golden_file)
        error("Golden file not found: $golden_file")
    end
    
    println("[LOAD] Loading MATLAB golden file: $golden_file")
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

function run_julia_simulation(Np, tmax)
    """
    Run Julia simulation with EXACT MATLAB parameters.
    
    MATLAB default parameters (from parse_main_args.m):
    - Kn = 1.0
    - Ma = 0.0
    - CFL = 0.5
    - rhol = 1.0
    - rhor = 0.01
    - T = 1.0
    - r110, r101, r011 = 0.0
    - flag2D = 0 (3D)
    """
    
    println("\n[RUN] Running Julia simulation with EXACT MATLAB parameters:")
    println("  Np = $Np")
    println("  tmax = $tmax")
    println("  Kn = 1.0, Ma = 0.0, CFL = 0.5 (MATLAB defaults)")
    println("  rhol = 1.0, rhor = 0.01 (MATLAB IC)")
    
    # Setup parameters EXACTLY as MATLAB does
    dx = 1.0 / Np
    dy = 1.0 / Np
    Nmom = 35
    nnmax = 20000000
    Kn = 1.0
    dtmax = Kn
    
    # Initial condition parameters - EXACTLY as MATLAB
    rhol = 1.0
    rhor = 0.01
    T = 1.0
    r110 = 0.0
    r101 = 0.0
    r011 = 0.0
    
    # Physical parameters - EXACTLY as MATLAB
    Ma = 0.0
    flag2D = 0
    CFL = 0.5
    
    # Diagnostic parameters
    symmetry_check_interval = 10
    
    # Package parameters
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
        enable_memory_tracking=false,
        debug_output=false
    )
    
    # Run simulation
    start_time = time()
    M_final, final_time, time_steps, grid_out = HyQMOM.simulation_runner(params)
    elapsed = time() - start_time
    
    println("  ✓ Julia simulation complete in $(round(elapsed, digits=2))s")
    println("  ✓ Julia final time: $final_time ($(time_steps) steps)")
    if M_final !== nothing
        println("  ✓ Moment array shape: ", size(M_final))
    end
    
    return (M=M_final, final_time=final_time, time_steps=time_steps, grid=grid_out)
end

function compare_results(matlab_result, julia_result; 
                        tolerance_abs=INTEGRATION_TOL_ABS, 
                        tolerance_rel=INTEGRATION_TOL_REL,
                        verbose=true)
    """Compare MATLAB and Julia results with detailed diagnostics."""
    
    M_matlab = matlab_result.M
    M_julia = julia_result.M
    
    if M_julia === nothing
        if verbose
            println("  ⚠ Julia result is nothing (not on rank 0)")
        end
        return false
    end
    
    if verbose
        println("\n" * "="^70)
        println("[COMPARE] COMPARISON RESULTS")
        println("="^70)
    end
    
    # Check dimensions
    if verbose
        println("\n1. Dimension Check:")
    end
    if size(M_matlab) != size(M_julia)
        if verbose
            println("  ✗ DIMENSION MISMATCH!")
            println("     MATLAB: ", size(M_matlab))
            println("     Julia:  ", size(M_julia))
        end
        return false
    elseif verbose
        println("  ✓ Dimensions match: ", size(M_julia))
    end
    
    # Check time
    if verbose
        println("\n2. Time Check:")
    end
    time_diff = abs(matlab_result.final_time - julia_result.final_time)
    if time_diff > 1e-10
        if verbose
            println("  ⚠ Time mismatch:")
            println("     MATLAB: $(matlab_result.final_time) ($(matlab_result.time_steps) steps)")
            println("     Julia:  $(julia_result.final_time) ($(julia_result.time_steps) steps)")
            println("     Diff:   $time_diff")
        end
    elseif verbose
        println("  ✓ Time matches: $(julia_result.final_time)")
        println("     MATLAB steps: $(matlab_result.time_steps)")
        println("     Julia steps:  $(julia_result.time_steps)")
    end
    
    # Check for NaN or Inf
    if verbose
        println("\n3. Numerical Health Check:")
    end
    nan_count = sum(isnan, M_julia)
    inf_count = sum(isinf, M_julia)
    
    if nan_count > 0 || inf_count > 0
        if verbose
            println("  ✗ NUMERICAL ISSUES:")
            println("     NaN count: $nan_count")
            println("     Inf count: $inf_count")
        end
        return false
    elseif verbose
        println("  ✓ No NaN or Inf values")
    end
    
    # Compute differences
    if verbose
        println("\n4. Value Comparison:")
    end
    diff = M_julia .- M_matlab
    abs_diff = abs.(diff)
    rel_diff = abs_diff ./ (abs.(M_matlab) .+ 1e-30)
    
    max_abs_diff = maximum(abs_diff)
    max_rel_diff = maximum(rel_diff)
    mean_abs_diff = sum(abs_diff) / length(abs_diff)
    mean_rel_diff = sum(rel_diff) / length(rel_diff)
    
    if verbose
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
        
        # Moment-by-moment analysis
        println("\n5. Moment-by-Moment Analysis:")
        println("-"^70)
        nx, ny, nmom = size(M_julia)
        
        moment_names = ["rho", "rhou", "rhov", "rhow", "rhoE", "M110", "M101", "M011", 
                        "M200", "M020", "M002"]
        
        for k in 1:min(11, nmom)
            moment_abs_diff = abs_diff[:, :, k]
            moment_rel_diff = rel_diff[:, :, k]
            
            max_abs = maximum(moment_abs_diff)
            max_rel = maximum(moment_rel_diff)
            mean_abs = sum(moment_abs_diff) / length(moment_abs_diff)
            
            moment_name = k <= length(moment_names) ? moment_names[k] : "M[$k]"
            
            @printf("  %-6s: max_abs=%.6e, max_rel=%.6e, mean_abs=%.6e", 
                    moment_name, max_abs, max_rel, mean_abs)
            
            if max_abs > tolerance_abs || max_rel > tolerance_rel
                print("  ⚠")
            end
            println()
        end
    end
    
    # Overall pass/fail
    passed = (max_abs_diff < tolerance_abs && max_rel_diff < tolerance_rel)
    
    if verbose
        println("\n" * "="^70)
        println("6. Final Verdict:")
        println("-"^70)
        @printf("  Tolerance (absolute): %.6e\n", tolerance_abs)
        @printf("  Tolerance (relative): %.6e\n", tolerance_rel)
        @printf("  Max absolute diff:    %.6e\n", max_abs_diff)
        @printf("  Max relative diff:    %.6e\n", max_rel_diff)
        println()
        
        if passed
            println("  ✓ TEST PASSED!")
            println("  Julia matches MATLAB within tolerance.")
        else
            println("  ✗ TEST FAILED!")
            println("  Differences exceed tolerance thresholds.")
            
            exceed_abs = sum(abs_diff .> tolerance_abs)
            exceed_rel = sum(rel_diff .> tolerance_rel)
            total_points = length(abs_diff)
            
            @printf("\n  Points exceeding absolute tolerance: %d / %d (%.2f%%)\n", 
                    exceed_abs, total_points, 100*exceed_abs/total_points)
            @printf("  Points exceeding relative tolerance: %d / %d (%.2f%%)\n",
                    exceed_rel, total_points, 100*exceed_rel/total_points)
        end
        
        println("="^70)
    end
    
    return passed
end

# Main test execution
function run_integration_tests()
    if !RUN_INTEGRATION
        @info "Skipping integration tests (TEST_INTEGRATION=false)"
        return
    end
    
    if !isfile(GOLDEN_FILE)
        @warn "MATLAB golden file not found: $GOLDEN_FILE"
        @warn "Run create_goldenfiles('ci') in MATLAB to generate it."
        @warn "Skipping integration tests."
        return
    end
    
    # Initialize MPI if needed
    mpi_state = ensure_mpi_initialized()
    
    if STANDALONE && mpi_state.rank == 0
        println("="^70)
        println("INTEGRATION TEST - Julia vs MATLAB")
        println("="^70)
        println("\nMPI Configuration: $(mpi_state.nprocs) rank(s)")
        println()
    end
    
    # Load MATLAB golden file (only rank 0)
    matlab_result = nothing
    if mpi_state.rank == 0
        try
            matlab_result = load_matlab_golden(GOLDEN_FILE)
        catch e
            @error "Error loading golden file" exception=e
            if HAS_MPI && MPI.Initialized()
                MPI.Finalize()
            end
            return
        end
    end
    
    # Run Julia simulation
    Np = 20
    tmax = 0.1
    
    julia_result = run_julia_simulation(Np, tmax)
    
    # Compare results (only on rank 0)
    if mpi_state.rank == 0 && matlab_result !== nothing
        test_passed = compare_results(matlab_result, julia_result, verbose=STANDALONE)
        
        if STANDALONE
            if HAS_MPI && MPI.Initialized()
                MPI.Finalize()
            end
            exit(test_passed ? 0 : 1)
        end
    end
    
    if HAS_MPI && MPI.Initialized() && STANDALONE
        MPI.Finalize()
    end
end

# If running within Test framework
if !STANDALONE
    @testset "Integration Tests" begin
        if !RUN_INTEGRATION
            @info "Skipping integration tests (TEST_INTEGRATION=false)"
        elseif !isfile(GOLDEN_FILE)
            @warn "MATLAB golden file not found: $GOLDEN_FILE"
            @test_skip "Integration test (no golden file)"
        else
            # Ensure MPI is initialized
            mpi_state = ensure_mpi_initialized()
            
            @testset "Load MATLAB Golden File" begin
                @test isfile(GOLDEN_FILE)
                
                global matlab_data = matread(GOLDEN_FILE)
                @test haskey(matlab_data, "golden_data")
                
                global golden_data = matlab_data["golden_data"]
                @test haskey(golden_data, "moments")
                @test haskey(golden_data, "parameters")
                
                global M_matlab = golden_data["moments"]["M"]
                global params_matlab = golden_data["parameters"]
                
                @test size(M_matlab, 1) == params_matlab["Np"]
                @test size(M_matlab, 2) == params_matlab["Np"]
                @test size(M_matlab, 3) == 35
            end
            
            @testset "Run Julia Simulation" begin
                Np = Int(params_matlab["Np"])
                tmax = params_matlab["tmax"]
                
                # Setup parameters
                params = (
                    Np=Np, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
                    dx=1.0/Np, dy=1.0/Np, Nmom=35, nnmax=20000000, dtmax=1.0,
                    rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
                    symmetry_check_interval=10,
                    enable_memory_tracking=false,
                    debug_output=false
                )
                
                global M_julia, final_time_julia, time_steps_julia, grid_julia = 
                    HyQMOM.simulation_runner(params)
                
                if M_julia !== nothing
                    @test size(M_julia) == size(M_matlab)
                    @test time_steps_julia == params_matlab["time_steps"]
                    @test abs(final_time_julia - params_matlab["final_time"]) < 1e-10
                end
            end
            
            @testset "Compare Results" begin
                if M_julia === nothing
                    @test_skip "Comparison (not on rank 0)"
                else
                    @test !any(isnan, M_julia)
                    @test !any(isinf, M_julia)
                    
                    # Compute differences
                    diff = M_julia .- M_matlab
                    abs_diff = abs.(diff)
                    rel_diff = abs_diff ./ (abs.(M_matlab) .+ 1e-30)
                    
                    max_abs_diff = maximum(abs_diff)
                    max_rel_diff = maximum(rel_diff)
                    
                    @test max_abs_diff < INTEGRATION_TOL_ABS
                    @test max_rel_diff < INTEGRATION_TOL_REL
                    
                    if max_abs_diff < INTEGRATION_TOL_ABS && max_rel_diff < INTEGRATION_TOL_REL
                        println("  ✓ Julia matches MATLAB within tolerance")
                    end
                end
            end
        end
    end
else
    # Running standalone
    run_integration_tests()
end

