"""
Find where the jet regions actually are in the IC
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
    tmax = 0.0,
    Kn = 1.0,
    Ma = 1.0,
    flag2D = 0,
    CFL = 0.7,
    dx = 1.0/40,
    dy = 1.0/40,
    dz = 1.0/20,
    Nmom = 35,
    nnmax = 0,
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

M_final, final_time, time_steps, grid = simulation_runner(params)

if rank == 0 && M_final !== nothing
    Np = params.Np
    Nz = params.Nz
    
    rho = M_final[:, :, :, 1]
    U = M_final[:, :, :, 2] ./ rho
    V = M_final[:, :, :, 6] ./ rho
    
    println("\n" * "="^70)
    println("FINDING JET REGIONS")
    println("="^70)
    
    # Find all points with high density (> 0.5)
    k_mid = div(Nz, 2)
    
    println("\nSearching middle z-plane (k=$k_mid) for high density (>0.5)...")
    
    jet_points = []
    for i in 1:Np
        for j in 1:Np
            if rho[i, j, k_mid] > 0.5
                push!(jet_points, (i, j, rho[i, j, k_mid], U[i, j, k_mid], V[i, j, k_mid]))
            end
        end
    end
    
    if length(jet_points) == 0
        println("\n❌ NO HIGH DENSITY REGIONS FOUND!")
        println("\nDensity statistics:")
        println(@sprintf("  Min: %.6f", minimum(rho)))
        println(@sprintf("  Max: %.6f", maximum(rho)))
        println(@sprintf("  Mean: %.6f", sum(rho)/length(rho)))
        
        # Show expected jet indices based on code
        Csize = floor(Int, 0.1 * Np)
        Mint = div(Np, 2) + 1
        Maxt = div(Np, 2) + 1 + Csize
        Minb = div(Np, 2) - Csize
        Maxb = div(Np, 2)
        
        println("\nExpected jet regions (from code):")
        println("  Csize = $Csize (10% of Np=$Np)")
        println("  Bottom jet: i ∈ [$Minb, $Maxb], j ∈ [$Minb, $Maxb]")
        println("  Top jet:    i ∈ [$Mint, $Maxt], j ∈ [$Mint, $Maxt]")
        
        println("\nChecking those specific regions:")
        i_test = div(Minb + Maxb, 2)
        j_test = div(Minb + Maxb, 2)
        println(@sprintf("  Bottom jet center (i=%d, j=%d): ρ=%.6f, U=%+.6f, V=%+.6f", 
                        i_test, j_test, rho[i_test, j_test, k_mid], 
                        U[i_test, j_test, k_mid], V[i_test, j_test, k_mid]))
        
        i_test = div(Mint + Maxt, 2)
        j_test = div(Mint + Maxt, 2)
        println(@sprintf("  Top jet center (i=%d, j=%d): ρ=%.6f, U=%+.6f, V=%+.6f", 
                        i_test, j_test, rho[i_test, j_test, k_mid], 
                        U[i_test, j_test, k_mid], V[i_test, j_test, k_mid]))
    else
        println("\n✓ Found $(length(jet_points)) high-density grid points")
        println("\nJet regions (first 10 points):")
        for (idx, (i, j, ρ, u, v)) in enumerate(jet_points[1:min(10, length(jet_points))])
            println(@sprintf("  Point %2d: (i=%2d, j=%2d) ρ=%.4f, U=%+.6f, V=%+.6f", 
                            idx, i, j, ρ, u, v))
        end
        
        # Analyze velocity patterns
        println("\nVelocity analysis:")
        positive_u = count(p -> p[4] > 0.5, jet_points)
        negative_u = count(p -> p[4] < -0.5, jet_points)
        positive_v = count(p -> p[5] > 0.5, jet_points)
        negative_v = count(p -> p[5] < -0.5, jet_points)
        
        println("  Points with U > +0.5: $positive_u")
        println("  Points with U < -0.5: $negative_u")
        println("  Points with V > +0.5: $positive_v")
        println("  Points with V < -0.5: $negative_v")
    end
    
    println("\n" * "="^70)
end

MPI.Finalize()

