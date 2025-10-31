"""
Configure MPI.jl to use system MPI binary.

Usage:
    module load openmpi/4.1.5  # or your MPI module
    julia --project=. scripts/setup_mpi.jl
"""

using Pkg

println("="^70)
println("MPI Configuration")
println("="^70)
println()

# Add MPIPreferences if not already present
if !haskey(Pkg.project().dependencies, "MPIPreferences")
    println("Adding MPIPreferences to Project.toml...")
    Pkg.add("MPIPreferences")
    println("  ✓ MPIPreferences added")
    println()
end

# Now use it to configure MPI
println("Configuring MPI.jl to use system binary...")
using MPIPreferences
MPIPreferences.use_system_binary()

println()
println("="^70)
println("✓ MPI configured for system binary")
println("="^70)
println()
println("LocalPreferences.toml contents:")
if isfile("LocalPreferences.toml")
    println(read("LocalPreferences.toml", String))
else
    println("  ERROR: LocalPreferences.toml not created!")
    println("  Make sure MPI module is loaded before running this script.")
end
println("="^70)

