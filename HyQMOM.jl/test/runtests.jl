"""
Main Test Entry Point for HyQMOM.jl

This file is the entry point for Julia's package testing system (Pkg.test()).
It runs all unit tests and optionally integration tests.

For more control over test execution, use the shell scripts:
- ./test/run_tests.sh         # All tests with nice formatting
- ./test/run_mpi_tests.sh     # MPI-specific tests

Environment Variables:
- TEST_INTEGRATION: Set to "false" to skip integration tests
"""

using Test
using HyQMOM
using LinearAlgebra

# Test tolerance
const TOL = 1e-10

@testset "HyQMOM.jl" begin
    @testset "Unit Tests" begin
        # Core functionality unit tests
        include("test_autogen.jl")
        include("test_moment_conversions.jl")
        include("test_initialization.jl")
        include("test_realizability.jl")
        include("test_closures.jl")
        include("test_numerical_schemes.jl")
        include("test_golden_files.jl")
    end
    
    # Integration tests (Julia vs MATLAB golden files)
    # The test file handles its own skip logic if golden files are missing
    # or if TEST_INTEGRATION=false
    include("test_integration.jl")
end
