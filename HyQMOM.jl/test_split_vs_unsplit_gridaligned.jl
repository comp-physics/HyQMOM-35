#!/usr/bin/env julia
"""
Test: Do split and unsplit methods give identical results for grid-aligned ICs?

For grid-aligned initial conditions (jets along X/Y/Z axes, no rotation),
both methods should theoretically give the same results since:
- Dimensional splitting: X/Y/Z flux updates
- 3D unsplit: Rotates to face normals, but for cardinal directions (±X,±Y,±Z),
  the rotation is a no-op (or simple swap), so it should reduce to the same fluxes.

This test checks if they match.
"""

using MPI
using HyQMOM
using Printf

MPI.Init()

rank = MPI.Comm_rank(MPI.COMM_WORLD)

# Small, fast test case with grid-aligned jets
base_params = (
    Nx = 16,
    Ny = 16,
    Nz = 16,
    tmax = 0.02,  # Short simulation
    CFL = 0.5,
    Kn = 1.0,
    Ma = 0.0,
    Pr = 1.0,
    flag2D = 0,
    
    xmin = -0.5,
    xmax = 0.5,
    ymin = -0.5,
    ymax = 0.5,
    zmin = -0.5,
    zmax = 0.5,
    dx = 1.0/16,
    dy = 1.0/16,
    dz = 1.0/16,
    
    Nmom = 35,
    nnmax = 10000,
    dtmax = 0.01,
    
    # Grid-aligned jets (default IC from simulation_runner)
    rhol = 1.0,
    rhor = 0.01,
    T = 1.0,
    r110 = 0.0,
    r101 = 0.0,
    r011 = 0.0,
    
    symmetry_check_interval = 1000,
    enable_memory_tracking = false,
    debug_output = false,
    homogeneous_x = false,
    homogeneous_y = false,
    homogeneous_z = false,
    
    output_interval = 0,
    snapshot_interval = 0
)

if rank == 0
    println("\n" * "="^70)
    println("TESTING: Split vs Unsplit for Grid-Aligned ICs")
    println("="^70)
    println("Question: Do they give identical results for non-rotated ICs?")
    println("Grid: $(base_params.Nx)×$(base_params.Ny)×$(base_params.Nz)")
    println("ICs: Grid-aligned jets (default, no rotation)")
    println("="^70)
end

# Run with dimensional splitting
if rank == 0
    println("\n1️⃣  Running DIMENSIONAL SPLITTING...")
end
params_split = merge(base_params, (use_3d_unsplit = false,))
t_start_split = time()
M_split, t_final_split, steps_split, _ = HyQMOM.simulation_runner(params_split)
walltime_split = time() - t_start_split

if rank == 0
    println(@sprintf("   ✓ Completed in %.2f s, %d steps", walltime_split, steps_split))
end

# Run with 3D unsplit
if rank == 0
    println("\n2️⃣  Running 3D UNSPLIT...")
end
params_unsplit = merge(base_params, (use_3d_unsplit = true,))
t_start_unsplit = time()
M_unsplit, t_final_unsplit, steps_unsplit, _ = HyQMOM.simulation_runner(params_unsplit)
walltime_unsplit = time() - t_start_unsplit

if rank == 0
    println(@sprintf("   ✓ Completed in %.2f s, %d steps", walltime_unsplit, steps_unsplit))
end

# Compare results
if rank == 0
    println("\n" * "="^70)
    println("COMPARISON RESULTS")
    println("="^70)
    
    # Check if they took the same number of steps
    if steps_split == steps_unsplit
        println(@sprintf("✓ Same number of steps: %d", steps_split))
    else
        println(@sprintf("✗ Different steps: split=%d, unsplit=%d", steps_split, steps_unsplit))
    end
    
    # Compute differences for all 35 moments
    println("\nMoment-by-Moment Comparison:")
    println("-"^70)
    
    max_errors = Float64[]
    mean_errors = Float64[]
    
    for m in 1:35
        diff = M_unsplit[:,:,:,m] .- M_split[:,:,:,m]
        max_err = maximum(abs.(diff))
        mean_err = sum(abs.(diff)) / length(diff)
        push!(max_errors, max_err)
        push!(mean_errors, mean_err)
        
        if m <= 10  # Print first 10 moments in detail
            @printf("  Moment %2d: max_err = %.6e, mean_err = %.6e\n", m, max_err, mean_err)
        end
    end
    
    println("  ... (moments 11-35)")
    
    # Overall statistics
    println("\n" * "="^70)
    println("OVERALL STATISTICS")
    println("="^70)
    
    overall_max = maximum(max_errors)
    overall_mean = sum(mean_errors) / length(mean_errors)
    
    @printf("Maximum error (all moments): %.6e\n", overall_max)
    @printf("Mean error (all moments):    %.6e\n", overall_mean)
    
    # Determine if they're "identical"
    tol_identical = 1e-12  # Machine precision tolerance
    tol_close = 1e-8       # Very close tolerance
    
    println("\n" * "="^70)
    println("VERDICT")
    println("="^70)
    
    if overall_max < tol_identical
        println("✅ IDENTICAL (machine precision)")
        println("   The methods give the same results for grid-aligned ICs!")
        println("   This confirms that 3D unsplit reduces to dimensional splitting")
        println("   for cardinal directions.")
    elseif overall_max < tol_close
        println("✅ NEARLY IDENTICAL (very close)")
        println("   Differences are extremely small (< 1e-8)")
        println("   Likely due to floating-point roundoff or order of operations.")
    else
        println("❌ DIFFERENT")
        println("   The methods give different results even for grid-aligned ICs!")
        println("   Max error: $(overall_max)")
        println("\n   Possible reasons:")
        println("   1. Implementation bug in 3D unsplit")
        println("   2. Order of operations differs (e.g., realizability timing)")
        println("   3. Flux computation differs even for cardinal directions")
    end
    
    # Performance comparison
    println("\n" * "="^70)
    println("PERFORMANCE")
    println("="^70)
    speedup = walltime_unsplit / walltime_split
    @printf("Dimensional Splitting: %.2f s\n", walltime_split)
    @printf("3D Unsplit:           %.2f s (%.2fx slower)\n", walltime_unsplit, speedup)
    
    println("="^70)
end

MPI.Finalize()

