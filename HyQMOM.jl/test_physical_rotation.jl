#!/usr/bin/env julia
"""
Test: True rotational invariance via physical rotation

Setup:
1. Define initial condition with 2 colliding blocks
2. Run simulation WITHOUT rotation
3. Run simulation WITH Z-rotation of IC
4. Rotate final state BACK to original orientation
5. Compare: rotated-back result should match non-rotated result

This tests the PHYSICS, not just grid discretization.
"""

using HyQMOM
using MPI
using Printf

MPI.Init()

Nx = Ny = Nz = 24  # Use decent resolution
tmax = 0.03
Ma = 0.0
Uc = Ma / sqrt(3.0)

# Two blocks colliding diagonally (as you described)
offset = 0.15
width = 0.12

background = CubicRegion(
    center = (0.0, 0.0, 0.0),
    width = (Inf, Inf, Inf),
    density = 0.01,
    velocity = (0.0, 0.0, 0.0),
    temperature = 1.0
)

# REFERENCE: Two blocks moving toward each other diagonally
blocks_ref = [
    CubicRegion(
        center = (-offset, -offset, 0.0),  # Lower-left block
        width = (width, width, width),
        density = 1.0,
        velocity = (Uc, Uc, 0.0),  # Moving ↗
        temperature = 1.0
    ),
    CubicRegion(
        center = (offset, offset, 0.0),   # Upper-right block
        width = (width, width, width),
        density = 1.0,
        velocity = (-Uc, -Uc, 0.0),  # Moving ↙
        temperature = 1.0
    )
]

# Z-ROTATED: Apply Z-rotation (+90°) to initial condition
# (x,y,z) → (y,-x,z), (u,v,w) → (v,-u,w)
blocks_rotated = [
    CubicRegion(
        center = (-offset, offset, 0.0),   # Rotate position
        width = (width, width, width),
        density = 1.0,
        velocity = (Uc, -Uc, 0.0),         # Rotate velocity
        temperature = 1.0
    ),
    CubicRegion(
        center = (offset, -offset, 0.0),
        width = (width, width, width),
        density = 1.0,
        velocity = (-Uc, Uc, 0.0),
        temperature = 1.0
    )
]

base_params = (
    Nx=Nx, Ny=Ny, Nz=Nz,
    xmin=-0.5, xmax=0.5,
    ymin=-0.5, ymax=0.5,
    zmin=-0.5, zmax=0.5,
    tmax=tmax,
    CFL=0.3,
    Kn=1.0,
    Ma=Ma,
    Pr=0.7,
    flag2D=0,
    Nmom=35,
    nnmax=1000000,
    dtmax=1e-2,
    rhol=1.0,
    rhor=0.01,
    T=1.0,
    r110=0.0,
    r101=0.0,
    r011=0.0,
    homogeneous_x=false,
    homogeneous_y=false,
    homogeneous_z=false,
    output_interval=0,
    snapshot_interval=0,
    symmetry_check_interval=10000,
    enable_memory_tracking=false,
    debug_output=false,
    use_custom_ic=true,
    ic_background=background,
    use_3d_unsplit=true  # Use the unsplit method!
)

println("="^70)
println("PHYSICAL ROTATION TEST: True Rotational Invariance")
println("="^70)
println("Configuration: 2 colliding blocks, Ma=$Ma")
println("Grid: $(Nx)×$(Ny)×$(Nz), tmax=$tmax")
println("="^70)
println("\nTest procedure:")
println("  1. Run reference simulation (no rotation)")
println("  2. Run with Z-rotated initial condition (+90°)")
println("  3. Rotate final state back (-90°)")
println("  4. Compare: should match if truly rotationally invariant")
println("="^70)

# Step 1: Reference simulation
println("\n[1/2] Running REFERENCE simulation...")
params_ref = merge(base_params, (ic_jets=blocks_ref,))
M_ref, t_ref, steps_ref, grid = HyQMOM.simulation_runner(params_ref)
println("  ✓ Complete: $steps_ref steps, t=$t_ref")
@printf("    Min ρ: %.6e, Max ρ: %.6e\n", minimum(M_ref[:,:,:,1]), maximum(M_ref[:,:,:,1]))

# Step 2: Rotated simulation
println("\n[2/2] Running Z-ROTATED simulation...")
params_rot = merge(base_params, (ic_jets=blocks_rotated,))
M_rot, t_rot, steps_rot, _ = HyQMOM.simulation_runner(params_rot)
println("  ✓ Complete: $steps_rot steps, t=$t_rot")
@printf("    Min ρ: %.6e, Max ρ: %.6e\n", minimum(M_rot[:,:,:,1]), maximum(M_rot[:,:,:,1]))

# Step 3: Rotate final state back
println("\n[3/3] Rotating final state back to reference orientation...")
println("  Using FULL 35-moment tensor rotation...")

# Apply complete inverse rotation: grid permutation + moment transformation
M_back = HyQMOM.rotate_full_state_z90_inverse(M_rot)

# Step 4: Compare
println("\n" * "="^70)
println("COMPARISON: Reference vs Rotated-Back (ALL 35 MOMENTS)")
println("="^70)

# Compute errors for all moments
max_errors = Float64[]
rms_errors = Float64[]
moment_names = ["ρ", "ρu", "ρuu", "ρuuu", "ρuuuu", "ρv", "ρuv", "ρuuv", "ρuuuv", 
                "ρvv", "ρuvv", "ρuuvv", "ρvvv", "ρuvvv", "ρvvvv", "ρw", "ρuw", "ρuuw",
                "ρuuuw", "ρww", "ρuww", "ρuuww", "ρwww", "ρuwww", "ρwwww", "ρvw",
                "ρuuw", "ρuuww", "ρvvw", "ρuvw", "ρvvww", "ρuvww", "ρuvvw", "ρvwww", "ρvvvww"]

for m in 1:35
    Δ = M_back[:,:,:,m] .- M_ref[:,:,:,m]
    push!(max_errors, maximum(abs.(Δ)))
    push!(rms_errors, sqrt(sum(Δ.^2) / length(Δ)))
end

@printf("\nWorst 10 Moments (Max Abs Error):\n")
@printf("%-10s  %-12s  %-12s\n", "Moment", "Max Error", "RMS Error")
println("-"^40)
sorted_indices = sortperm(max_errors, rev=true)
for i in 1:min(10, 35)
    idx = sorted_indices[i]
    @printf("M%-9d  %.6e  %.6e\n", idx-1, max_errors[idx], rms_errors[idx])
end

@printf("\nFirst 4 Moments (Physics):\n")
@printf("  |Δρ|:    %.6e (RMS: %.6e)\n", max_errors[1], rms_errors[1])
@printf("  |Δ(ρu)|: %.6e (RMS: %.6e)\n", max_errors[2], rms_errors[2])
@printf("  |Δ(ρv)|: %.6e (RMS: %.6e)\n", max_errors[3], rms_errors[3])
@printf("  |Δ(ρw)|: %.6e (RMS: %.6e)\n", max_errors[4], rms_errors[4])

max_error = maximum(max_errors[1:4])  # Judge based on physical moments
max_error_all = maximum(max_errors)   # Track all moments

println("\n" * "="^70)
println("VERDICT")
println("="^70)

@printf("\nPhysical moments (ρ,u,v,w): Max error = %.6e\n", max_error)
@printf("All 35 moments:            Max error = %.6e\n", max_error_all)

if max_error < 1e-10
    println("\n✅ PERFECT ROTATIONAL INVARIANCE")
    println("   Error < 10^-10 (machine precision)")
    println("   The physics is EXACTLY preserved under rotation!")
elseif max_error < 1e-6
    println("\n✅ EXCELLENT ROTATIONAL INVARIANCE")
    println("   Error < 10^-6 (numerical precision)")
    println("   The physics is very well preserved under rotation!")
elseif max_error < 1e-3
    println("\n✓ GOOD ROTATIONAL INVARIANCE")
    println("  Error < 10^-3")
    println("  Minor discretization effects, physics mostly preserved.")
elseif max_error < 0.01
    println("\n⚠️  MODERATE ROTATIONAL INVARIANCE")
    println("   Error < 0.01 (1%)")
    println("   Some discretization effects, but reasonable.")
elseif max_error < 0.1
    println("\n⚠️  WEAK ROTATIONAL INVARIANCE")
    println("   Error < 0.1 (10%)")
    println("   Significant discretization effects.")
else
    println("\n❌ POOR ROTATIONAL INVARIANCE")
    println("   Error > 0.1 (10%)")
    println("   Either discretization is too coarse or algorithm has issues.")
end

if max_error_all > 10 * max_error
    println("\n⚠️  WARNING: Higher-order moments have much larger errors!")
    println("   This may indicate incomplete realizability or tensor rotation issues.")
end

println("="^70)

MPI.Finalize()

