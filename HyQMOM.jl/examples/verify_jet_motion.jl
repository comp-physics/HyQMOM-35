"""
Verify that jets are moving in the correct directions

This runs a very short simulation and checks that:
1. Initially: bottom-left has (+U,+V), top-right has (-U,-V)
2. After time evolution: jets have moved toward each other
"""

using HyQMOM
using MPI
using Printf

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

params = (
    Np = 40,
    Nz = 20,
    tmax = 0.01,  # Very short time
    Kn = 1.0,
    Ma = 1.0,
    flag2D = 0,
    CFL = 0.5,
    dx = 1.0/40,
    dy = 1.0/40,
    dz = 1.0/20,
    Nmom = 35,
    nnmax = 5,  # Just a few steps
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
    println("Running short simulation to verify jet motion...")
    println("Jets should move TOWARD each other (crossing pattern)")
end

M_final, final_time, time_steps, grid = simulation_runner(params)

if rank == 0 && M_final !== nothing
    Np = params.Np
    Nz = params.Nz
    k_mid = div(Nz, 2)
    
    rho = M_final[:, :, :, 1]
    U = M_final[:, :, :, 2] ./ rho
    V = M_final[:, :, :, 6] ./ rho
    
    println("\n" * "="^70)
    println("JET MOTION VERIFICATION (after $time_steps steps)")
    println("="^70)
    
    # Find bottom-left jet (should have +U, +V and be moving toward top-right)
    # Look in the left-bottom quadrant
    function find_max_density_location(rho_slice, i_range, j_range)
        max_rho = 0.0
        best_i, best_j = i_range[1], j_range[1]
        for i in i_range
            for j in j_range
                if rho_slice[i, j] > max_rho
                    max_rho = rho_slice[i, j]
                    best_i, best_j = i, j
                end
            end
        end
        return best_i, best_j, max_rho
    end
    
    i_bl, j_bl, max_rho_bl = find_max_density_location(rho[:, :, k_mid], 1:div(Np,2), 1:div(Np,2))
    
    println("\nBottom-Left Jet (highest density in lower-left quadrant):")
    println(@sprintf("  Position: (i=%d, j=%d)", i_bl, j_bl))
    println(@sprintf("  Density: %.4f", rho[i_bl, j_bl, k_mid]))
    println(@sprintf("  U velocity: %+.6f", U[i_bl, j_bl, k_mid]))
    println(@sprintf("  V velocity: %+.6f", V[i_bl, j_bl, k_mid]))
    
    if U[i_bl, j_bl, k_mid] > 0.5 && V[i_bl, j_bl, k_mid] > 0.5
        println("  ✓ CORRECT: Moving toward TOP-RIGHT (+U, +V)")
    else
        println("  ✗ WRONG: Velocities don't indicate top-right motion!")
    end
    
    # Find top-right jet (should have -U, -V and be moving toward bottom-left)
    i_tr, j_tr, max_rho_tr = find_max_density_location(rho[:, :, k_mid], div(Np,2)+1:Np, div(Np,2)+1:Np)
    
    println("\nTop-Right Jet (highest density in upper-right quadrant):")
    println(@sprintf("  Position: (i=%d, j=%d)", i_tr, j_tr))
    println(@sprintf("  Density: %.4f", rho[i_tr, j_tr, k_mid]))
    println(@sprintf("  U velocity: %+.6f", U[i_tr, j_tr, k_mid]))
    println(@sprintf("  V velocity: %+.6f", V[i_tr, j_tr, k_mid]))
    
    if U[i_tr, j_tr, k_mid] < -0.5 && V[i_tr, j_tr, k_mid] < -0.5
        println("  ✓ CORRECT: Moving toward BOTTOM-LEFT (-U, -V)")
    else
        println("  ✗ WRONG: Velocities don't indicate bottom-left motion!")
    end
    
    # Check if jets have moved closer
    println("\n" * "="^70)
    println("MOTION CHECK")
    println("="^70)
    println(@sprintf("Bottom-left jet at (i=%d, j=%d)", i_bl, j_bl))
    println(@sprintf("Top-right jet at   (i=%d, j=%d)", i_tr, j_tr))
    
    # Initial positions should be around i~17, j~17 for BL and i~23, j~23 for TR
    # After motion toward center, BL should move right/up, TR should move left/down
    
    if i_bl > 17 || j_bl > 17
        println("✓ Bottom-left jet has moved RIGHT and/or UP (toward center)")
    else
        println("⚠  Bottom-left jet position seems unchanged or moved away")
    end
    
    if i_tr < 23 || j_tr < 23
        println("✓ Top-right jet has moved LEFT and/or DOWN (toward center)")
    else
        println("⚠  Top-right jet position seems unchanged or moved away")
    end
    
    println("\n" * "="^70)
    println("CONCLUSION")
    println("="^70)
    println("Expected behavior:")
    println("  • Bottom-left jet: (+U, +V) moving ↗ toward top-right")
    println("  • Top-right jet:   (-U, -V) moving ↙ toward bottom-left")
    println("  • Result: Jets CROSS in the center")
    println("="^70)
end

MPI.Finalize()

