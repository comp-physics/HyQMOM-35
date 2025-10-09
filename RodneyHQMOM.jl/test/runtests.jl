using Test
using RodneyHQMOM
using LinearAlgebra

# Test tolerance
const TOL = 1e-10

# Golden file configuration (must be at module level, not in testset)
const GOLDEN_FILE = joinpath(@__DIR__, "..", "..", "goldenfiles", 
                              "goldenfile_mpi_1ranks_Np20_tmax100.mat")
const RUN_GOLDEN_TEST = get(ENV, "TEST_GOLDEN_SIMULATION", "true") != "false"

@testset "RodneyHQMOM.jl" begin
    # Unit tests
    include("test_autogen.jl")
    include("test_moment_conversions.jl")
    include("test_initialization.jl")
    include("test_realizability.jl")
    include("test_closures.jl")
    include("test_numerical_schemes.jl")
    include("test_golden_files.jl")
    include("test_integration_1rank.jl")
    
    # Full simulation golden file test (Julia vs MATLAB)
    # Only run if golden file exists and TEST_GOLDEN_SIMULATION is not set to "false"
    if RUN_GOLDEN_TEST && isfile(GOLDEN_FILE)
        @testset "MATLAB Golden File (Full Simulation)" begin
            println("\n" * "="^70)
            println("Testing Julia vs MATLAB Golden File")
            println("="^70)
            
            # Run the golden file test
            # This is a standalone test that includes its own MPI initialization
            # For Pkg.test(), we run it without MPI
            include("test_matlab_golden_simple.jl")
        end
    elseif RUN_GOLDEN_TEST && !isfile(GOLDEN_FILE)
        @warn "MATLAB golden file not found: $GOLDEN_FILE\nSkipping full simulation test."
    else
        @info "Skipping MATLAB golden file test (TEST_GOLDEN_SIMULATION=false)"
    end
end
