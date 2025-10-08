using Test
using RodneyHQMOM
using LinearAlgebra

# Test tolerance
const TOL = 1e-10

@testset "RodneyHQMOM.jl" begin
    include("test_autogen.jl")
    include("test_moment_conversions.jl")
    include("test_initialization.jl")
    include("test_realizability.jl")
    include("test_closures.jl")
    include("test_numerical_schemes.jl")
    include("test_golden_files.jl")
    include("test_integration_1rank.jl")
end
