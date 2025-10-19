#!/usr/bin/env julia
"""
Test: Do split and unsplit give same results for crossing jets (the working case)?
"""

using MPI
using HyQMOM
using Printf

MPI.Init()

rank = MPI.Comm_rank(MPI.COMM_WORLD)

# Use the SAME configuration as run_3d_custom_jets.jl crossing case
Ma = 0.0
Kn = 1.0
rhol = 1.0
rhor = 0.01
T = 1.0

Lx = Ly = Lz = 1.0
jet_width = 0.1 * Lx
Uc = Ma / sqrt(3.0)
offset = jet_width * 0.6

background = HyQMOM.CubicRegion(
    center = (0.0, 0.0, 0.0),
    width = (Inf, Inf, Inf),
    density = rhor,
    velocity = (0.0, 0.0, 0.0),
    temperature = T
)

jets = [
    HyQMOM.CubicRegion(
        center = (-offset, -offset, -offset),
        width = (jet_width, jet_width, jet_width),
        density = rhol,
        velocity = (Uc, Uc, Uc),
        temperature = T
    ),
    HyQMOM.CubicRegion(
        center = (offset, offset, offset),
        width = (jet_width, jet_width, jet_width),
        density = rhol,
        velocity = (-Uc, -Uc, -Uc),
        temperature = T
    ),
]

# Base parameters
base_params = (
    Nx = 16,
    Ny = 16,
    Nz = 16,
    tmax = 0.01,  # Very short test
    CFL = 0.5,
    Kn = Kn,
    Ma = Ma,
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
    
    rhol = rhol,
    rhor = rhor,
    T = T,
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
    snapshot_interval = 0,
    
    # Use custom IC!
    use_custom_ic = true,
    ic_background = background,
    ic_jets = jets
)

if rank == 0
    println("\n" * "="^70)
    println("TESTING: Split vs Unsplit for Crossing Jets")
    println("="^70)
    println("Using the SAME config as run_3d_custom_jets.jl")
    println("Grid: $(base_params.Nx)×$(base_params.Ny)×$(base_params.Nz)")
    println("="^70)
end

# Test 1: Dimensional Splitting
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

# Test 2: 3D Unsplit
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

# Compare
if rank == 0
    println("\n" * "="^70)
    println("COMPARISON")
    println("="^70)
    
    max_diff_rho = maximum(abs.(M_unsplit[:,:,:,1] .- M_split[:,:,:,1]))
    @printf("Max |Δρ|: %.6e\n", max_diff_rho)
    
    if max_diff_rho < 1e-10
        println("✅ IDENTICAL!")
    elseif max_diff_rho < 1e-6
        println("✅ VERY CLOSE")
    else
        println("❌ DIFFERENT")
    end
    
    println("="^70)
end

MPI.Finalize()

