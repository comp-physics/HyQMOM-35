#!/usr/bin/env julia
"""
Test 3D unsplit method (rotationally invariant algorithm).

This script runs a simple 3D simulation using the new 3D unsplit method,
which maintains rotational invariance by rotating moment vectors to align
with face normals before flux computation.

Usage:
    mpiexec -n 1 julia --project=. examples/test_3d_unsplit.jl
    
Note: This method is ~2x slower than dimensional splitting but provides
      perfect rotational invariance.
"""

using MPI
using HyQMOM

MPI.Init()

# Simulation parameters
params = (
    Nx = 20,
    Ny = 20,
    Nz = 20,
    tmax = 0.1,
    CFL = 0.5,
    Kn = 1.0,
    Ma = 0.0,
    Pr = 1.0,
    flag2D = 0,
    
    # Grid spacing
    xmin = -0.5,
    xmax = 0.5,
    ymin = -0.5,
    ymax = 0.5,
    zmin = -0.5,
    zmax = 0.5,
    dx = 1.0/20,
    dy = 1.0/20,
    dz = 1.0/20,
    
    # Moment parameters
    Nmom = 35,
    nnmax = 10000,
    dtmax = 0.01,
    
    # IC parameters
    rhol = 1.0,
    rhor = 0.01,
    T = 1.0,
    r110 = 0.0,
    r101 = 0.0,
    r011 = 0.0,
    
    # Diagnostics
    symmetry_check_interval = 10,
    enable_memory_tracking = false,
    debug_output = false,
    homogeneous_x = false,
    homogeneous_y = false,
    homogeneous_z = false,
    
    # Output
    output_interval = 0,
    snapshot_interval = 0,
    
    # Algorithm selection: USE 3D UNSPLIT (rotationally invariant)
    use_3d_unsplit = true
)

rank = MPI.Comm_rank(MPI.COMM_WORLD)

if rank == 0
    println("="^70)
    println("3D UNSPLIT METHOD TEST")
    println("="^70)
    println("Grid: $(params.Nx)×$(params.Ny)×$(params.Nz)")
    println("Ma=$(params.Ma), Kn=$(params.Kn), tmax=$(params.tmax)")
    println("Method: 3D Unsplit (Rotationally Invariant)")
    println("="^70)
end

# Run simulation
M_final, t_final, steps, grid = HyQMOM.simulation_runner(params)

if rank == 0
    println("\n" * "="^70)
    println("SIMULATION COMPLETE")
    println("="^70)
    println("Final time: $t_final")
    println("Time steps: $steps")
    println("Final density range: [$(minimum(M_final[:,:,:,1])), $(maximum(M_final[:,:,:,1]))]")
    println("="^70)
end

MPI.Finalize()

