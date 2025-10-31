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

# Add MPIPreferences if not already present
if !haskey(Pkg.project().dependencies, "MPIPreferences")
    println("Adding MPIPreferences...")
    Pkg.add("MPIPreferences")
    println()
end

using MPIPreferences

println("Configuring MPI.jl to use system binary...")
MPIPreferences.use_system_binary()

println()
println("✓ MPI configured for system binary")
println()
println("LocalPreferences.toml:")
if isfile("LocalPreferences.toml")
    println(read("LocalPreferences.toml", String))
else
    println("  (file not created - may indicate an error)")
end
println("="^70)

