#!/usr/bin/env julia
"""
Compare dimensional splitting vs 3D unsplit methods.

This script runs both methods and compares:
1. Final state differences
2. Performance (walltime)
3. Conservation properties

Usage:
    mpiexec -n 1 julia --project=. examples/compare_methods.jl
"""

using MPI
using HyQMOM
using Printf

MPI.Init()

# Base parameters
base_params = (
    Nx = 16,
    Ny = 16,
    Nz = 16,
    tmax = 0.05,
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

rank = MPI.Comm_rank(MPI.COMM_WORLD)

if rank == 0
    println("\n" * "="^70)
    println("METHOD COMPARISON: Dimensional Splitting vs 3D Unsplit")
    println("="^70)
    println("Grid: $(base_params.Nx)×$(base_params.Ny)×$(base_params.Nz)")
    println("tmax: $(base_params.tmax), Ma: $(base_params.Ma), Kn: $(base_params.Kn)")
end

# Test 1: Dimensional Splitting
if rank == 0
    println("\n" * "-"^70)
    println("TEST 1: Dimensional Splitting (MATLAB-verified)")
    println("-"^70)
end

params_split = merge(base_params, (use_3d_unsplit = false,))
t_start_split = time()
M_split, t_final_split, steps_split, _ = HyQMOM.simulation_runner(params_split)
walltime_split = time() - t_start_split

if rank == 0
    println(@sprintf("  Completed: %.2f s walltime, %d steps", walltime_split, steps_split))
    println(@sprintf("  ρ range: [%.6e, %.6e]", minimum(M_split[:,:,:,1]), maximum(M_split[:,:,:,1])))
end

# Test 2: 3D Unsplit
if rank == 0
    println("\n" * "-"^70)
    println("TEST 2: 3D Unsplit (Rotationally Invariant)")
    println("-"^70)
end

params_unsplit = merge(base_params, (use_3d_unsplit = true,))
t_start_unsplit = time()
M_unsplit, t_final_unsplit, steps_unsplit, _ = HyQMOM.simulation_runner(params_unsplit)
walltime_unsplit = time() - t_start_unsplit

if rank == 0
    println(@sprintf("  Completed: %.2f s walltime, %d steps", walltime_unsplit, steps_unsplit))
    println(@sprintf("  ρ range: [%.6e, %.6e]", minimum(M_unsplit[:,:,:,1]), maximum(M_unsplit[:,:,:,1])))
end

# Comparison
if rank == 0
    println("\n" * "="^70)
    println("COMPARISON SUMMARY")
    println("="^70)
    
    # Performance
    speedup = walltime_unsplit / walltime_split
    println(@sprintf("Performance:"))
    println(@sprintf("  Splitting:  %.2f s", walltime_split))
    println(@sprintf("  Unsplit:    %.2f s (%.2fx slower)", walltime_unsplit, speedup))
    
    # Differences
    diff_rho = M_unsplit[:,:,:,1] .- M_split[:,:,:,1]
    diff_u = M_unsplit[:,:,:,2] .- M_split[:,:,:,2]
    diff_v = M_unsplit[:,:,:,3] .- M_split[:,:,:,3]
    diff_w = M_unsplit[:,:,:,4] .- M_split[:,:,:,4]
    
    println(@sprintf("\nState Differences:"))
    println(@sprintf("  Max |Δρ|: %.6e", maximum(abs.(diff_rho))))
    println(@sprintf("  Max |Δu|: %.6e", maximum(abs.(diff_u))))
    println(@sprintf("  Max |Δv|: %.6e", maximum(abs.(diff_v))))
    println(@sprintf("  Max |Δw|: %.6e", maximum(abs.(diff_w))))
    
    # Mass conservation
    mass_split = sum(M_split[:,:,:,1])
    mass_unsplit = sum(M_unsplit[:,:,:,1])
    println(@sprintf("\nMass Conservation:"))
    println(@sprintf("  Splitting: %.12e", mass_split))
    println(@sprintf("  Unsplit:   %.12e", mass_unsplit))
    println(@sprintf("  Difference: %.6e", abs(mass_unsplit - mass_split)))
    
    println("="^70)
    println("\nNote: Methods are expected to give different results for most ICs")
    println("      due to dimensional splitting anisotropy. The 3D unsplit method")
    println("      maintains rotational invariance at the cost of ~2x performance.")
    println("="^70)
end

MPI.Finalize()

