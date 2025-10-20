#!/usr/bin/env julia
"""
Validation test for rotational invariance of dimensional splitting method.

This test quantifies the rotational anisotropy by:
1. Running a reference simulation (no rotation)
2. Running with a Z-rotated initial condition (+90°)
3. Rotating the final state back (-90°)
4. Comparing to measure anisotropy

Expected result: < 0.2% error for physical moments (ρ, ρu, ρv, ρw)
"""

using HyQMOM
using MPI
using Printf

MPI.Init()

rank = MPI.Comm_rank(MPI.COMM_WORLD)
nprocs = MPI.Comm_size(MPI.COMM_WORLD)

# Test parameters - adjust as needed
Nx = 40  # Use 40 for thorough test, 20 for quick test
Ny = 40
Nz = 40
tmax = 0.1
Ma = 0.0
Uc = 0.0  # Stationary for Ma=0

if rank == 0
    println("="^70)
    println("ROTATIONAL INVARIANCE TEST (Dimensional Splitting)")
    println("="^70)
    println("Grid: $(Nx)×$(Ny)×$(Nz)")
    println("Processors: $nprocs")
    println("tmax: $tmax")
    println("Ma: $Ma")
    println("="^70)
end

background = CubicRegion(center=(0.0, 0.0, 0.0), width=(Inf, Inf, Inf),
                        density=0.01, velocity=(0.0, 0.0, 0.0), temperature=1.0)

# Reference: Two jets colliding diagonally
offset = 0.15
jets_ref = [
    CubicRegion(center=(-offset, -offset, 0.0), width=(0.12, 0.12, 0.12),
                density=1.0, velocity=(Uc, Uc, 0.0), temperature=1.0),
    CubicRegion(center=(offset, offset, 0.0), width=(0.12, 0.12, 0.12),
                density=1.0, velocity=(-Uc, -Uc, 0.0), temperature=1.0)
]

# Z-rotated (+90°): (x,y) → (y,-x), (u,v) → (v,-u)
jets_rot = [
    CubicRegion(center=(-offset, offset, 0.0), width=(0.12, 0.12, 0.12),
                density=1.0, velocity=(Uc, -Uc, 0.0), temperature=1.0),
    CubicRegion(center=(offset, -offset, 0.0), width=(0.12, 0.12, 0.12),
                density=1.0, velocity=(-Uc, Uc, 0.0), temperature=1.0)
]

base = (Nx=Nx, Ny=Ny, Nz=Nz, xmin=-0.5, xmax=0.5, ymin=-0.5, ymax=0.5, zmin=-0.5, zmax=0.5,
        tmax=tmax, CFL=0.3, Kn=1.0, Ma=Ma, Pr=0.7, flag2D=0, Nmom=35, nnmax=1000000,
        dtmax=1e-2, rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
        homogeneous_x=false, homogeneous_y=false, homogeneous_z=false,
        output_interval=0, snapshot_interval=0, symmetry_check_interval=10000,
        enable_memory_tracking=false, debug_output=false, use_custom_ic=true,
        ic_background=background)

if rank == 0
    println("\n[1/2] Running REFERENCE simulation...")
end
params_ref = merge(base, (ic_jets=jets_ref,))
t0 = time()
M_ref, t_ref, steps_ref, grid = HyQMOM.simulation_runner(params_ref)
t_ref_wall = time() - t0

if rank == 0
    println("  ✓ Complete: $steps_ref steps, t=$t_ref, wall time=$(round(t_ref_wall,digits=1))s")
    @printf("    Min ρ: %.6e, Max ρ: %.6e\n", minimum(M_ref[:,:,:,1]), maximum(M_ref[:,:,:,1]))
end

if rank == 0
    println("\n[2/2] Running Z-ROTATED simulation...")
end
params_rot = merge(base, (ic_jets=jets_rot,))
t0 = time()
M_rot, t_rot, steps_rot, _ = HyQMOM.simulation_runner(params_rot)
t_rot_wall = time() - t0

if rank == 0
    println("  ✓ Complete: $steps_rot steps, t=$t_rot, wall time=$(round(t_rot_wall,digits=1))s")
    @printf("    Min ρ: %.6e, Max ρ: %.6e\n", minimum(M_rot[:,:,:,1]), maximum(M_rot[:,:,:,1]))
end

if rank == 0
    println("\n[3/3] Rotating final state back and comparing...")
    
    # Simple Z-90° rotation (inverse): (x,y) → (-y,x), (u,v) → (-v,u)
    # Grid permutation: M_back[i,j,k,:] = M_rot[j,Ny+1-i,k,:]
    M_back = zeros(size(M_rot))
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        j_rot = i
        i_rot = Ny + 1 - j
        # Moments: [ρ, ρu, ρv, ρw, ...]
        M_back[i,j,k,1] = M_rot[i_rot,j_rot,k,1]   # ρ unchanged
        M_back[i,j,k,2] = -M_rot[i_rot,j_rot,k,3]  # ρu = -ρv_rot
        M_back[i,j,k,3] = M_rot[i_rot,j_rot,k,2]   # ρv = ρu_rot
        M_back[i,j,k,4] = M_rot[i_rot,j_rot,k,4]   # ρw unchanged
        # Copy remaining moments (not rotationally invariant, but for mass check)
        for m in 5:35
            M_back[i,j,k,m] = M_rot[i_rot,j_rot,k,m]
        end
    end
    
    # Compute errors for physical moments
    Δρ = M_back[:,:,:,1] .- M_ref[:,:,:,1]
    Δρu = M_back[:,:,:,2] .- M_ref[:,:,:,2]
    Δρv = M_back[:,:,:,3] .- M_ref[:,:,:,3]
    Δρw = M_back[:,:,:,4] .- M_ref[:,:,:,4]
    
    max_err_ρ = maximum(abs.(Δρ))
    max_err_u = maximum(abs.(Δρu))
    max_err_v = maximum(abs.(Δρv))
    max_err_w = maximum(abs.(Δρw))
    max_err = max(max_err_ρ, max_err_u, max_err_v, max_err_w)
    
    # Mass conservation check
    mass_ref = sum(M_ref[:,:,:,1]) * (grid.dx * grid.dy * grid.dz)
    mass_rot = sum(M_rot[:,:,:,1]) * (grid.dx * grid.dy * grid.dz)
    mass_back = sum(M_back[:,:,:,1]) * (grid.dx * grid.dy * grid.dz)
    
    println("\n" * "="^70)
    println("RESULTS")
    println("="^70)
    
    @printf("\nPhysical Moments (Rotational Anisotropy):\n")
    @printf("  |Δρ|:    %.6e (%.2f%%)\n", max_err_ρ, max_err_ρ*100)
    @printf("  |Δ(ρu)|: %.6e (%.2f%%)\n", max_err_u, max_err_u*100)
    @printf("  |Δ(ρv)|: %.6e (%.2f%%)\n", max_err_v, max_err_v*100)
    @printf("  |Δ(ρw)|: %.6e (%.2f%%)\n", max_err_w, max_err_w*100)
    @printf("\nMaximum error: %.6e (%.2f%%)\n", max_err, max_err*100)
    
    @printf("\nMass Conservation:\n")
    @printf("  Reference:     %.12e\n", mass_ref)
    @printf("  Rotated:       %.12e\n", mass_rot)
    @printf("  Rotated-back:  %.12e\n", mass_back)
    @printf("  Ref vs Rot:    %.6e (%.4f%%)\n", abs(mass_rot - mass_ref), abs(mass_rot - mass_ref)/mass_ref*100)
    @printf("  Ref vs Back:   %.6e (%.4f%%)\n", abs(mass_back - mass_ref), abs(mass_back - mass_ref)/mass_ref*100)
    
    println("\n" * "="^70)
    
    if max_err < 0.002
        println("✅ PASS: < 0.2% rotational anisotropy")
        println("   This is excellent for a Cartesian grid method!")
    elseif max_err < 0.01
        println("✅ PASS: < 1% rotational anisotropy")
        println("   This is very good for a Cartesian grid method")
    elseif max_err < 0.05
        println("⚠️  ACCEPTABLE: < 5% rotational anisotropy")
        println("   Within typical range for dimensional splitting")
    else
        println("❌ FAIL: > 5% rotational anisotropy")
        println("   Investigate potential issues")
    end
    
    println("="^70)
end

MPI.Finalize()

