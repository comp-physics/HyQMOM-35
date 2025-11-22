"""
CI-friendly test runner for Julia 1.9

This script works around Julia 1.9's Test extension precompilation issues.
Use this in CI instead of `Pkg.test()` directly.

Usage:
    julia --project=. test/ci_test.jl
"""

# Workaround for Julia 1.9 Test extension issues
# Must be set before any package loading
if VERSION < v"1.10"
    ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
end

using Pkg

# Ensure we're in the project directory
cd(dirname(@__DIR__))

# Activate the test environment
Pkg.activate(".")

# Instantiate without auto-precompiling extensions
println("Instantiating packages (Julia $(VERSION))...")
Pkg.instantiate()

# Now run the actual tests
println("\nRunning tests...")
include("runtests.jl")

