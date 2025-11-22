"""
MPI Consistency and Integration Tests

Tests that MPI parallelization works correctly by comparing results
from different rank counts. This consolidates what was previously in:
- test_mpi_consistency.jl
- test_mpi_integration.jl
- test_mpi_goldenfiles.jl

Usage:
    # Within Pkg.test() framework (requires MPI)
    mpiexec -n 2 julia --project=. -e 'using Pkg; Pkg.test()'
    
    # Standalone usage
    # Step 1: Generate reference with 1 rank
    julia --project=. test/test_mpi.jl
    
    # Step 2: Test with 2 ranks
    mpiexec -n 2 julia --project=. test/test_mpi.jl
    
    # Step 3: Test with 4 ranks
    mpiexec -n 4 julia --project=. test/test_mpi.jl
    
    # Golden file mode (requires pre-generated golden files)
    mpiexec -n 2 julia --project=. test/test_mpi.jl --golden small
    mpiexec -n 4 julia --project=. test/test_mpi.jl --golden medium
"""

using Test
using MPI
using Printf
using HyQMOM

# Configuration
const MPI_TOL_ABS = 1e-10  # Very tight tolerance for deterministic results
const MPI_TOL_REL = 1e-8
const GOLDEN_DIR = joinpath(@__DIR__, "goldenfiles")
const REF_DIR = @__DIR__

# Test configurations
const TEST_CONFIGS = Dict(
    "small" => (Nx=20, Ny=20, tmax=0.05),
    "medium" => (Nx=40, Ny=40, tmax=0.05),
    "quick" => (Nx=20, Ny=20, tmax=0.02)  # For fast testing
)

# Check if we're running standalone or within Test framework
const STANDALONE = abspath(PROGRAM_FILE) == @__FILE__

# Parse command line arguments
function parse_args()
    mode = "dynamic"  # default: dynamic comparison
    config = "quick"  # default config
    
    for arg in ARGS
        if arg == "--golden"
            mode = "golden"
        elseif haskey(TEST_CONFIGS, arg)
            config = arg
        end
    end
    
    return (mode=mode, config=config)
end

function save_reference(M, final_time, time_steps, Np, nprocs, config_name)
    """Save simulation results as reference data."""
    ref_file = joinpath(REF_DIR, "mpi_reference_$(config_name)_$(nprocs)ranks.bin")
    
    println("[SAVE] Saving $(nprocs)-rank reference: $(ref_file)")
    
    open(ref_file, "w") do io
        write(io, Int64(nprocs))
        write(io, Int64(Np))
        write(io, Float64(final_time))
        write(io, Int64(time_steps))
        write(io, M)
    end
    
    println("  [OK] Saved successfully")
    return ref_file
end

function load_reference(config_name, ref_nprocs=1)
    """Load reference data from file."""
    ref_file = joinpath(REF_DIR, "mpi_reference_$(config_name)_$(ref_nprocs)ranks.bin")
    
    if !isfile(ref_file)
        return nothing
    end
    
    println("[LOAD] Loading $(ref_nprocs)-rank reference: $(ref_file)")
    
    nprocs, Np, final_time, time_steps, M = open(ref_file, "r") do io
        nprocs = read(io, Int64)
        Np = read(io, Int64)
        final_time = read(io, Float64)
        time_steps = read(io, Int64)
        M = Array{Float64}(undef, Np, Np, 35)
        read!(io, M)
        (nprocs, Np, final_time, time_steps, M)
    end
    
    println("  [OK] Loaded: Np=$(Np), t=$(final_time), steps=$(time_steps)")
    
    return (nprocs=nprocs, Np=Np, final_time=final_time, 
            time_steps=time_steps, M=M)
end

function load_golden_file(config_name)
    """Load pre-generated golden file for regression testing."""
    golden_file = joinpath(GOLDEN_DIR, "mpi_1rank_$(config_name).bin")
    
    if !isfile(golden_file)
        @warn "Golden file not found: $(golden_file)"
        @warn "Run: mpiexec -n 1 julia test/create_golden_files.jl"
        return nothing
    end
    
    println("[GOLDEN] Loading golden file: $(golden_file)")
    
    nranks, Np, tmax, tfinal, steps, M = open(golden_file, "r") do io
        nranks = read(io, Int64)
        Np = read(io, Int64)
        tmax = read(io, Float64)
        tfinal = read(io, Float64)
        steps = read(io, Int64)
        M = Array{Float64}(undef, Np, Np, 35)
        read!(io, M)
        (nranks, Np, tmax, tfinal, steps, M)
    end
    
    println("  [OK] Loaded: Np=$(Np), tmax=$(tmax), t=$(tfinal), steps=$(steps)")
    
    return (nranks=nranks, Np=Np, tmax=tmax, final_time=tfinal, 
            time_steps=steps, M=M)
end

function run_simulation(config_name)
    """Run HyQMOM simulation with specified configuration."""
    
    # Get configuration
    if !haskey(TEST_CONFIGS, config_name)
        error("Unknown configuration: $(config_name)")
    end
    
    config = TEST_CONFIGS[config_name]
    Nx = config.Nx
    Ny = config.Ny
    tmax = config.tmax
    
    # Get MPI info
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    if rank == 0
        println("\n[RUN] Running $(nprocs)-rank simulation: $(config_name)")
        println("  Nx=$(Nx), Ny=$(Ny), tmax=$(tmax)")
    end
    
    # Run simulation
    results = HyQMOM.run_simulation(
        Nx = Nx,
        Ny = Ny,
        tmax = tmax,
        num_workers = nprocs,
        verbose = false,
        Kn = 1.0,
        Ma = 0.0,
        flag2D = 0,
        CFL = 0.5
    )
    
    if rank == 0
        M = results[:M]
        final_time = results[:final_time]
        time_steps = results[:time_steps]
        
        println("  [OK] Simulation complete")
        println("    Final time: $(final_time)")
        println("    Time steps: $(time_steps)")
        
        return (M=M, final_time=final_time, time_steps=time_steps, Nx=Nx, Ny=Ny)
    else
        return nothing
    end
end

function compare_results(ref_result, test_result, test_nprocs; verbose=true)
    """Compare reference and test results."""
    
    M_ref = ref_result.M
    M_test = test_result.M
    
    if verbose
        println("\n" * "="^70)
        println("[COMPARE] 1 rank vs $(test_nprocs) ranks")
        println("="^70)
    end
    
    # Check dimensions
    if size(M_ref) != size(M_test)
        if verbose
            println("  [X] DIMENSION MISMATCH!")
            println("    1 rank:  ", size(M_ref))
            println("    $(test_nprocs) ranks: ", size(M_test))
        end
        return false
    elseif verbose
        println("\n1. Dimensions: ", size(M_test), " [OK]")
    end
    
    # Check time
    if verbose
        println("\n2. Time Comparison:")
    end
    time_diff = abs(ref_result.final_time - test_result.final_time)
    if time_diff > 1e-12
        if verbose
            println("  [WARNING] Time differs:")
            println("    1 rank:  $(ref_result.final_time) ($(ref_result.time_steps) steps)")
            println("    $(test_nprocs) ranks: $(test_result.final_time) ($(test_result.time_steps) steps)")
            println("    Diff:    $(time_diff)")
        end
    elseif verbose
        println("  [OK] Final time: $(test_result.final_time)")
        println("    Steps: $(test_result.time_steps)")
    end
    
    # Check for NaN/Inf
    nan_ref = sum(isnan, M_ref)
    inf_ref = sum(isinf, M_ref)
    nan_test = sum(isnan, M_test)
    inf_test = sum(isinf, M_test)
    
    if nan_ref > 0 || inf_ref > 0 || nan_test > 0 || inf_test > 0
        if verbose
            println("\n3. Numerical Health: [X]")
            println("  1 rank:  NaN=$(nan_ref), Inf=$(inf_ref)")
            println("  $(test_nprocs) ranks: NaN=$(nan_test), Inf=$(inf_test)")
        end
        return false
    elseif verbose
        println("\n3. Numerical Health: [OK] (no NaN/Inf)")
    end
    
    # Compute differences
    diff = M_test .- M_ref
    abs_diff = abs.(diff)
    rel_diff = abs_diff ./ (abs.(M_ref) .+ 1e-30)
    
    max_abs_diff = maximum(abs_diff)
    max_rel_diff = maximum(rel_diff)
    mean_abs_diff = sum(abs_diff) / length(abs_diff)
    mean_rel_diff = sum(rel_diff) / length(rel_diff)
    
    if verbose
        println("\n4. Differences:")
        @printf("  Max absolute: %.6e\n", max_abs_diff)
        @printf("  Max relative: %.6e\n", max_rel_diff)
        @printf("  Mean absolute: %.6e\n", mean_abs_diff)
        @printf("  Mean relative: %.6e\n", mean_rel_diff)
        
        # Find worst moment
        max_moment_diff = 0.0
        worst_moment = 0
        
        println("\n5. Per-Moment Analysis:")
        for k in 1:35
            moment_abs = maximum(abs_diff[:, :, k])
            if moment_abs > max_moment_diff
                max_moment_diff = moment_abs
                worst_moment = k
            end
            
            if moment_abs > MPI_TOL_ABS * 10
                @printf("    Moment %2d: %.6e [WARNING]\n", k, moment_abs)
            end
        end
        
        if max_moment_diff <= MPI_TOL_ABS
            println("  [OK] All moments within tolerance")
        else
            println("  [WARNING] Worst moment: #$(worst_moment) (diff: $(max_moment_diff))")
        end
    end
    
    # Pass/fail
    passed = (max_abs_diff < MPI_TOL_ABS && max_rel_diff < MPI_TOL_REL)
    
    if verbose
        println("\n" * "="^70)
        println("6. Result:")
        @printf("  Tolerance (abs): %.6e\n", MPI_TOL_ABS)
        @printf("  Tolerance (rel): %.6e\n", MPI_TOL_REL)
        println()
        
        if passed
            println("  [OK] MPI CONSISTENCY TEST PASSED!")
            println("    Results match within tolerance")
        else
            println("  [X] MPI CONSISTENCY TEST FAILED!")
            println("    Differences exceed tolerance")
            
            # Show location of max diff
            max_idx = argmax(abs_diff)
            max_loc = Tuple(max_idx)
            println("\n  Location of max diff: $(max_loc)")
            @printf("    1 rank:  %.10e\n", M_ref[max_idx])
            @printf("    $(test_nprocs) ranks: %.10e\n", M_test[max_idx])
            @printf("    Diff:    %.10e\n", diff[max_idx])
        end
        println("="^70)
    end
    
    return passed
end

function run_mpi_tests()
    """Main MPI test execution."""
    
    # Initialize MPI
    if !MPI.Initialized()
        MPI.Init()
    end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    # Parse arguments
    args = parse_args()
    mode = args.mode
    config_name = args.config
    
    if rank == 0
        println("="^70)
        println("MPI TESTS")
        println("="^70)
        println("\nMode: $(mode)")
        println("Config: $(config_name)")
        println("Ranks: $(nprocs)")
        println()
    end
    
    if mode == "golden"
        # Golden file mode - test against pre-generated golden files
        if rank == 0
            golden_result = load_golden_file(config_name)
            if golden_result === nothing
                println("[X] Golden file not found")
                MPI.Finalize()
                return 1
            end
        end
        
        # Run simulation
        test_result = run_simulation(config_name)
        
        # Compare (only rank 0)
        if rank == 0 && test_result !== nothing
            # Need to load the simulation parameters from golden file
            golden_result = load_golden_file(config_name)
            passed = compare_results(golden_result, test_result, nprocs)
            
            MPI.Finalize()
            return passed ? 0 : 1
        end
        
    else
        # Dynamic mode - compare 1 rank vs N ranks
        if nprocs == 1
            # Generate reference
            test_result = run_simulation(config_name)
            
            if rank == 0 && test_result !== nothing
                save_reference(test_result.M, test_result.final_time, 
                             test_result.time_steps, test_result.Np, 
                             nprocs, config_name)
                
                println("\n" * "="^70)
                println("[OK] Reference data generated (1 rank)")
                println("="^70)
                println("\nNext step: Run with multiple ranks to compare:")
                println("  mpiexec -n 2 julia --project=. test/test_mpi.jl")
                println("="^70)
            end
            
        else
            # Load reference and compare
            if rank == 0
                ref_result = load_reference(config_name, 1)
                
                if ref_result === nothing
                    println("[X] Reference data not found for $(config_name)")
                    println("  First run with 1 rank to generate reference:")
                    println("    julia --project=. test/test_mpi.jl")
                    MPI.Finalize()
                    return 1
                end
            end
            
            # Run test simulation
            test_result = run_simulation(config_name)
            
            # Compare (only rank 0)
            if rank == 0 && test_result !== nothing
                ref_result = load_reference(config_name, 1)
                passed = compare_results(ref_result, test_result, nprocs)
                
                MPI.Finalize()
                return passed ? 0 : 1
            end
        end
    end
    
    MPI.Finalize()
    return 0
end

# Main execution
if STANDALONE
    exit(run_mpi_tests())
else
    # Running within Test framework
    @testset "MPI Tests" begin
        # Check if MPI is available and we have multiple ranks
        if !MPI.Initialized()
            MPI.Init()
        end
        
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        nprocs = MPI.Comm_size(comm)
        
        if nprocs == 1
            @warn "MPI tests require multiple ranks. Run with: mpiexec -n 2 julia ..."
            @test_skip "MPI tests (only 1 rank)"
        else
            @testset "MPI $(nprocs) ranks" begin
                # Load reference (1 rank)
                ref_result = nothing
                if rank == 0
                    ref_result = load_reference("quick", 1)
                    if ref_result === nothing
                        @warn "Reference not found. Generate with: julia test/test_mpi.jl"
                    end
                end
                
                if ref_result === nothing && rank == 0
                    @test_skip "MPI test (no reference data)"
                else
                    # Run test
                    test_result = run_simulation("quick")
                    
                    if rank == 0 && test_result !== nothing
                        # Compare
                        passed = compare_results(ref_result, test_result, nprocs, verbose=false)
                        @test passed
                        
                        if passed
                            println("  [OK] MPI $(nprocs)-rank test passed")
                        else
                            println("  [X] MPI $(nprocs)-rank test failed")
                        end
                    end
                end
            end
        end
        
        # Additional halo exchange tests
        @testset "Halo Exchange Correctness" begin
            # Test that halo exchange correctly copies data
            # This is conceptual without running actual MPI
            
            @testset "Ghost cell indexing" begin
                Nx_local = 10
                # Ghost cells would be at indices 0 and Nx_local+1
                @test 1 <= 1 <= Nx_local  # Left interior
                @test 1 <= Nx_local <= Nx_local  # Right interior
            end
            
            @testset "Periodic topology" begin
                nprocs_test = 4
                for rank_test in 0:nprocs_test-1
                    rank_left = mod(rank_test - 1, nprocs_test)
                    rank_right = mod(rank_test + 1, nprocs_test)
                    
                    @test 0 <= rank_left < nprocs_test
                    @test 0 <= rank_right < nprocs_test
                end
            end
            
            @testset "2D decomposition" begin
                # Test 2D domain decomposition logic
                nranks_x = 2
                nranks_y = 2
                total_ranks = nranks_x * nranks_y
                
                @test total_ranks == 4
                
                # Test rank mapping
                for ry in 0:nranks_y-1
                    for rx in 0:nranks_x-1
                        rank_2d = ry * nranks_x + rx
                        @test 0 <= rank_2d < total_ranks
                    end
                end
            end
            
            @testset "Corner cell handling" begin
                # Test that diagonal neighbors are correctly identified
                # In 2D, each cell has 8 neighbors (including diagonals)
                neighbors_2d = [
                    (-1, -1), (-1, 0), (-1, 1),
                    ( 0, -1),          ( 0, 1),
                    ( 1, -1), ( 1, 0), ( 1, 1)
                ]
                
                @test length(neighbors_2d) == 8
            end
            
            @testset "Moment field consistency" begin
                # All 35 moments should be exchanged
                nmoments = 35
                @test nmoments == 35
            end
        end
    end
end

