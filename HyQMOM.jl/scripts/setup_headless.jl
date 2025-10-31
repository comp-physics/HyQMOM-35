"""
Setup script for headless HPC environments.
Run this BEFORE Pkg.instantiate() to remove visualization dependencies.

Usage:
    julia --project=. scripts/setup_headless.jl
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
"""

using Pkg

println("="^70)
println("HyQMOM Headless Setup")
println("="^70)
println("Removing visualization packages incompatible with headless systems...")
println()

# Packages that require OpenGL/X11
viz_packages = ["GLMakie", "FileIO", "ColorSchemes", "LaTeXStrings"]

for pkg in viz_packages
    try
        Pkg.rm(pkg; mode=PKGMODE_PROJECT)
        println("  ✓ Removed: $pkg")
    catch e
        println("  ⊘ Skipped: $pkg (not a direct dependency)")
    end
end

println()
println("="^70)
println("✓ Headless setup complete!")
println("="^70)
println()
println("Next steps:")
println("  1. julia --project=. -e 'using Pkg; Pkg.instantiate()'")
println("  2. Set: export HYQMOM_SKIP_PLOTTING=true")
println("  3. Configure MPI: julia --project=. scripts/setup_mpi.jl")
println()

