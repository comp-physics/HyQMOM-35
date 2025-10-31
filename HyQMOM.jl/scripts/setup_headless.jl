"""
Setup script for headless HPC environments.
Run this BEFORE Pkg.instantiate() to remove visualization dependencies.

Usage:
    julia --project=. scripts/setup_headless.jl
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
"""

using Pkg

# Disable auto-precompile to avoid stale cache errors during package removal
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"

println("="^70)
println("HyQMOM Headless Setup")
println("="^70)
println("Removing packages incompatible with headless systems...")
println()

# Packages to remove for headless compatibility:
# - GLMakie, FileIO, ColorSchemes, LaTeXStrings: require OpenGL/X11
# - MAT: pulls HDF5_jll with MPIExt → OpenMPI_jll conflicts with system MPI
remove_packages = [
    "GLMakie",
    "FileIO", 
    "ColorSchemes",
    "LaTeXStrings",
    "MAT",
    "HDF5"  # If present as explicit dep
]

println("Removing packages...")
for pkg in remove_packages
    try
        Pkg.rm(pkg; mode=PKGMODE_PROJECT)
        println("  ✓ Removed: $pkg")
    catch e
        println("  ⊘ Skipped: $pkg (not in dependencies)")
    end
end
println()

# Clear precompiled cache AFTER removing packages
# This avoids stale references to removed packages
println("Clearing precompiled cache...")
try
    depot_path = first(DEPOT_PATH)
    compiled_dir = joinpath(depot_path, "compiled", "v$(VERSION.major).$(VERSION.minor)")
    if isdir(compiled_dir)
        rm(compiled_dir, recursive=true, force=true)
        println("  ✓ Cache cleared")
    else
        println("  ⊘ No cache to clear")
    end
catch e
    println("  ⊘ Cache clear failed: $(typeof(e))")
end

# Note: We do NOT call Pkg.resolve() here because it can trigger precompilation
# with stale caches, causing "Dates not found" errors. The next Pkg.instantiate()
# will handle dependency resolution correctly.

println()
println("="^70)
println("✓ Headless setup complete!")
println("="^70)
println()
println("Removed packages:")
println("  - GLMakie, FileIO, ColorSchemes, LaTeXStrings (visualization)")
println("  - MAT, HDF5 (MPI JLL conflicts)")
println()
println("Next steps:")
println("  1. julia --project=. -e 'using Pkg; Pkg.instantiate()'")
println("  2. julia --project=. scripts/setup_mpi.jl")
println()
println("Note: JLD2 remains available for snapshot I/O")
println()

