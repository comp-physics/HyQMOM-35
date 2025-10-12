"""
Check initial condition velocities for crossing jets

This script verifies that:
- Bottom-left jet has velocity (+U, +V) moving toward top-right
- Top-right jet has velocity (-U, -V) moving toward bottom-left
"""

using HyQMOM
using MPI
using Printf

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

# Small grid for checking IC
params = (
    Np = 40,
    Nz = 20,
    tmax = 0.0,  # NO TIME EVOLUTION - just check IC
    Kn = 1.0,
    Ma = 1.0,
    flag2D = 0,
    CFL = 0.7,
    dx = 1.0/40,
    dy = 1.0/40,
    dz = 1.0/20,
    Nmom = 35,
    nnmax = 0,  # NO TIME STEPS
    dtmax = 1e-2,
    rhol = 1.0,
    rhor = 0.01,
    T = 1.0,
    r110 = 0.0,
    r101 = 0.0,
    r011 = 0.0,
    symmetry_check_interval = 100,
    homogeneous_z = false,
    debug_output = false,
    enable_memory_tracking = false
)

if rank == 0
    println("Checking initial condition velocities...")
    println("Ma = $(params.Ma), Uc = $(params.Ma / sqrt(2.0))")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

if rank == 0 && M_final !== nothing
    Np = params.Np
    Nz = params.Nz
    
    # Extract density and velocities
    rho = M_final[:, :, :, 1]
    U = M_final[:, :, :, 2] ./ rho
    V = M_final[:, :, :, 6] ./ rho
    W = M_final[:, :, :, 16] ./ rho
    
    println("\n" * "="^70)
    println("INITIAL CONDITION VERIFICATION")
    println("="^70)
    
    # Check bottom-left jet region (should have high density, +U, +V)
    i_bl = div(Np, 4)  # Quarter way in x
    j_bl = div(Np, 4)  # Quarter way in y
    k_mid = div(Nz, 2) # Middle in z
    
    println("\nBottom-Left Jet Region (i=$i_bl, j=$j_bl, k=$k_mid):")
    println(@sprintf("  Density: %.4f", rho[i_bl, j_bl, k_mid]))
    println(@sprintf("  U velocity: %+.4f (should be POSITIVE, ~%.4f)", U[i_bl, j_bl, k_mid], params.Ma/sqrt(2.0)))
    println(@sprintf("  V velocity: %+.4f (should be POSITIVE, ~%.4f)", V[i_bl, j_bl, k_mid], params.Ma/sqrt(2.0)))
    println(@sprintf("  W velocity: %+.4f (should be ~0)", W[i_bl, j_bl, k_mid]))
    
    if rho[i_bl, j_bl, k_mid] > 0.5
        if U[i_bl, j_bl, k_mid] > 0 && V[i_bl, j_bl, k_mid] > 0
            println("  ✓ CORRECT: Bottom-left jet moving toward TOP-RIGHT")
        else
            println("  ✗ WRONG: Velocities have incorrect signs!")
        end
    else
        println("  ⚠  WARNING: Low density, might not be in jet region")
    end
    
    # Check top-right jet region (should have high density, -U, -V)
    i_tr = div(3*Np, 4)  # Three-quarters way in x
    j_tr = div(3*Np, 4)  # Three-quarters way in y
    
    println("\nTop-Right Jet Region (i=$i_tr, j=$j_tr, k=$k_mid):")
    println(@sprintf("  Density: %.4f", rho[i_tr, j_tr, k_mid]))
    println(@sprintf("  U velocity: %+.4f (should be NEGATIVE, ~%.4f)", U[i_tr, j_tr, k_mid], -params.Ma/sqrt(2.0)))
    println(@sprintf("  V velocity: %+.4f (should be NEGATIVE, ~%.4f)", V[i_tr, j_tr, k_mid], -params.Ma/sqrt(2.0)))
    println(@sprintf("  W velocity: %+.4f (should be ~0)", W[i_tr, j_tr, k_mid]))
    
    if rho[i_tr, j_tr, k_mid] > 0.5
        if U[i_tr, j_tr, k_mid] < 0 && V[i_tr, j_tr, k_mid] < 0
            println("  ✓ CORRECT: Top-right jet moving toward BOTTOM-LEFT")
        else
            println("  ✗ WRONG: Velocities have incorrect signs!")
        end
    else
        println("  ⚠  WARNING: Low density, might not be in jet region")
    end
    
    # Check background region
    i_bg = div(Np, 2)
    j_bg = 1  # Edge
    
    println("\nBackground Region (i=$i_bg, j=$j_bg, k=$k_mid):")
    println(@sprintf("  Density: %.4f (should be ~%.2f)", rho[i_bg, j_bg, k_mid], params.rhor))
    println(@sprintf("  U velocity: %+.4f (should be ~0)", U[i_bg, j_bg, k_mid]))
    println(@sprintf("  V velocity: %+.4f (should be ~0)", V[i_bg, j_bg, k_mid]))
    
    println("\n" * "="^70)
    println("SUMMARY")
    println("="^70)
    println("Expected configuration:")
    println("  Bottom-left cube: ρ=1.0, U=+0.707, V=+0.707 → moving ↗ (toward top-right)")
    println("  Top-right cube:   ρ=1.0, U=-0.707, V=-0.707 → moving ↙ (toward bottom-left)")
    println("  Background:       ρ=0.01, U=0, V=0")
    println("\nThe jets should CROSS in the center of the domain.")
    println("="^70)
end

MPI.Finalize()

