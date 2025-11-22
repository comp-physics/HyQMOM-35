#!/usr/bin/env julia
# Restore full desktop dependencies to Project.toml
# Use this if you accidentally ran setup_headless.jl on your desktop,
# or want to restore visualization packages on a cluster.

println()
println("="^70)
println("Restoring Desktop Dependencies")
println("="^70)
println()

using Pkg

# Packages to restore
restore_packages = [
    "GLMakie",
    "FileIO", 
    "ColorSchemes",
    "LaTeXStrings"
]

println("Adding visualization packages...")
for pkg in restore_packages
    try
        Pkg.add(pkg)
        println("  ✓ Added: $pkg")
    catch e
        println("  ✗ Failed to add $pkg: ", e)
    end
end

println()
println("="^70)
println("✓ Desktop dependencies restored!")
println("="^70)
println()
println("Added packages:")
println("  - GLMakie, FileIO, ColorSchemes, LaTeXStrings (visualization)")
println()
println("Next steps:")
println("  1. julia --project=. -e 'using Pkg; Pkg.instantiate()'")
println("  2. julia --project=. -e 'using HyQMOM; println(\"✓ Ready!\")'")
println()
println("Note: MAT and HDF5 are NOT restored by default (cause MPI conflicts).")
println("      Add them manually if needed: Pkg.add([\"MAT\", \"HDF5\"])")
println()

