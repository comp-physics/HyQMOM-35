#!/usr/bin/env julia
# Load MATLAB step 4 state and compare with Julia step 5 computation

using RodneyHQMOM
using MAT
using Printf

println("="^70)
println("STEP 4→5 DETAILED COMPARISON: MATLAB vs JULIA")
println("="^70)

# Load MATLAB state at end of step 4
println("\n📥 Loading MATLAB step 4 state...")
matlab_data = matread("../matlab_step4_complete.mat")

M_matlab = matlab_data["M"]
t_matlab = matlab_data["t"]
dt_matlab = matlab_data["dt"]
nx = Int(matlab_data["nx"])
ny = Int(matlab_data["ny"])
halo = Int(matlab_data["halo"])
Nmom = Int(matlab_data["Nmom"])

println("  ✓ Loaded: nx=$nx, ny=$ny, halo=$halo, t=$t_matlab, dt=$dt_matlab")

# Parameters
Np = 20
Kn = 1.0
Ma = 0.0
flag2D = 0
CFL = 0.5
dx = 2.0 / Np
dy = 2.0 / Np

# Now run Julia from this state for ONE MORE STEP (step 5)
println("\n🔧 Setting up Julia simulation from MATLAB state...")
M = copy(M_matlab)  # Start from MATLAB's state

# Track cell (7,13)
track_i = 7
track_j = 13
ih_track = track_i + halo
jh_track = track_j + halo

println("\n" * "="^70)
println("STARTING STEP 5 FROM MATLAB STATE")
println("="^70)

@printf("\n📍 Cell (%d,%d) at START of step 5:\n", track_i, track_j)
@printf("  M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", 
        M[ih_track, jh_track, 1:5]...)

# Allocate arrays
Mnp = similar(M)
Fx = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
Fy = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)

vpxmin = zeros(Float64, nx, ny)
vpxmax = zeros(Float64, nx, ny)
vpymin = zeros(Float64, nx, ny)
vpymax = zeros(Float64, nx, ny)

println("\n" * "-"*70)
println("PHASE 1: Flux Computation")
println("-"*70)

# Compute fluxes for interior cells
for i in 1:nx
    for j in 1:ny
        ih = i + halo
        jh = j + halo
        MOM = M[ih, jh, :]
        
        _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
        vpxmin[i,j], vpxmax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
        vpymin[i,j], vpymax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
        Mx, My, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
        
        Fx[ih, jh, :] = Mx
        Fy[ih, jh, :] = My
        Mnp[ih, jh, :] = Mr
        
        # Detailed logging for tracked cell
        if i == track_i && j == track_j
            @printf("\nCell (%d,%d) after flux computation:\n", i, j)
            @printf("  Input MOM[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", MOM[1:5]...)
            @printf("  Output Mr[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mr[1:5]...)
            @printf("  Fx[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mx[1:5]...)
            @printf("  Fy[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", My[1:5]...)
        end
    end
end

M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]

@printf("\nCell (%d,%d) after bulk assignment:\n", track_i, track_j)
@printf("  M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", 
        M[ih_track, jh_track, 1:5]...)

# Apply boundary conditions (copy)
println("\nApplying boundary conditions...")
M[1:halo, :, :] = M[halo+1:halo+halo, :, :]
M[halo+nx+1:end, :, :] = M[halo+nx-halo+1:halo+nx, :, :]
M[:, 1:halo, :] = M[:, halo+1:halo+halo, :]
M[:, halo+ny+1:end, :] = M[:, halo+ny-halo+1:halo+ny, :]

# Exchange Fx, Fy boundaries
Fx[1:halo, :, :] = Fx[halo+1:halo+halo, :, :]
Fx[halo+nx+1:end, :, :] = Fx[halo+nx-halo+1:halo+nx, :, :]
Fx[:, 1:halo, :] = Fx[:, halo+1:halo+halo, :]
Fx[:, halo+ny+1:end, :] = Fx[:, halo+ny-halo+1:halo+ny, :]

Fy[1:halo, :, :] = Fy[halo+1:halo+halo, :, :]
Fy[halo+nx+1:end, :, :] = Fy[halo+nx-halo+1:halo+nx, :, :]
Fy[:, 1:halo, :] = Fy[:, halo+1:halo+halo, :]
Fy[:, halo+ny+1:end, :] = Fy[:, halo+ny-halo+1:halo+ny, :]

println("\n" * "-"*70)
println("PHASE 2: Time Step Computation")
println("-"*70)

vmax = maximum([abs.(vpxmin); abs.(vpxmax); abs.(vpymin); abs.(vpymax)])
dt = CFL * min(dx, dy) / vmax

@printf("vmax = %.6e\n", vmax)
@printf("dt = %.6e\n", dt)
@printf("MATLAB dt was: %.6e\n", dt_matlab)
@printf("Difference: %.6e\n", abs(dt - dt_matlab))

println("\n" * "-"*70)
println("PHASE 3: Flux Updates")
println("-"*70)

# Apply flux updates
for i in 1:nx
    for j in 1:ny
        ih = i + halo
        jh = j + halo
        
        # X-direction
        Mnpx = M[ih, jh, :] - dt/dx * (Fx[ih+1, jh, :] - Fx[ih, jh, :])
        
        # Y-direction  
        Mnpy = M[ih, jh, :] - dt/dy * (Fy[ih, jh+1, :] - Fy[ih, jh, :])
        
        # Combine
        M_combined = Mnpx + Mnpy - M[ih, jh, :]
        
        Mnp[ih, jh, :] = M_combined
        
        if i == track_i && j == track_j
            @printf("\nCell (%d,%d) flux update:\n", i, j)
            @printf("  Mnpx[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mnpx[1:5]...)
            @printf("  Mnpy[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mnpy[1:5]...)
            @printf("  Combined[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", M_combined[1:5]...)
        end
    end
end

M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]

@printf("\nCell (%d,%d) after flux update:\n", track_i, track_j)
@printf("  M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", 
        M[ih_track, jh_track, 1:5]...)

# Apply boundary conditions
M[1:halo, :, :] = M[halo+1:halo+halo, :, :]
M[halo+nx+1:end, :, :] = M[halo+nx-halo+1:halo+nx, :, :]
M[:, 1:halo, :] = M[:, halo+1:halo+halo, :]
M[:, halo+ny+1:end, :] = M[:, halo+ny-halo+1:halo+ny, :]

println("\n" * "-"*70)
println("PHASE 4: Realizability Enforcement") 
println("-"*70)

# Realizability
for i in 1:nx
    for j in 1:ny
        ih = i + halo
        jh = j + halo
        MOM = M[ih, jh, :]
        
        if i == track_i && j == track_j
            @printf("\n🔍 Cell (%d,%d) BEFORE realizability:\n", i, j)
            @printf("  Input[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", MOM[1:5]...)
        end
        
        _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma, debug_label="[realiz1]")
        
        if i == track_i && j == track_j
            @printf("  After call 1[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mr[1:5]...)
        end
        
        vpxmin_temp, vpxmax_temp, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
        
        if i == track_i && j == track_j
            @printf("  After eigen-x[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mr[1:5]...)
        end
        
        vpymin_temp, vpymax_temp, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
        
        if i == track_i && j == track_j
            @printf("  After eigen-y[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mr[1:5]...)
        end
        
        _, _, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma, debug_label="[realiz2]")
        
        if i == track_i && j == track_j
            @printf("  After call 2[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", Mr[1:5]...)
            
            # Check if exploded
            if abs(Mr[3]) > 1e6
                println("\n❌ ❌ ❌  EXPLOSION DETECTED IN REALIZABILITY! ❌ ❌ ❌")
                @printf("  Mr[3] = %.6e (way too large!)\n", Mr[3])
            end
        end
        
        Mnp[ih, jh, :] = Mr
    end
end

M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]

@printf("\n📍 Cell (%d,%d) after realizability:\n", track_i, track_j)
@printf("  M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", 
        M[ih_track, jh_track, 1:5]...)

println("\n" * "-"*70)
println("PHASE 5: Collision Operator")
println("-"*70)

# Collision
for i in 1:nx
    for j in 1:ny
        ih = i + halo
        jh = j + halo
        MOM = M[ih, jh, :]
        MMC = collision35(MOM, dt, Kn)
        Mnp[ih, jh, :] = MMC
        
        if i == track_i && j == track_j
            @printf("\nCell (%d,%d) collision:\n", i, j)
            @printf("  Before[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", MOM[1:5]...)
            @printf("  After[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", MMC[1:5]...)
        end
    end
end

M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]

# Final boundary conditions
M[1:halo, :, :] = M[halo+1:halo+halo, :, :]
M[halo+nx+1:end, :, :] = M[halo+nx-halo+1:halo+nx, :, :]
M[:, 1:halo, :] = M[:, halo+1:halo+halo, :]
M[:, halo+ny+1:end, :] = M[:, halo+ny-halo+1:halo+ny, :]

println("\n" * "="^70)
println("STEP 5 COMPLETE")
println("="^70)

@printf("\n📍 Cell (%d,%d) at END of step 5:\n", track_i, track_j)
@printf("  M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", 
        M[ih_track, jh_track, 1:5]...)

# Check for explosions
M_cell_final = M[ih_track, jh_track, :]
if any(abs.(M_cell_final) .> 1e6)
    println("\n❌ ❌ ❌  CELL EXPLODED! ❌ ❌ ❌")
    for (idx, val) in enumerate(M_cell_final)
        if abs(val) > 1e6
            @printf("  M[%d] = %.6e\n", idx, val)
        end
    end
elseif any(isnan.(M_cell_final))
    println("\n❌ ❌ ❌  NaN DETECTED! ❌ ❌ ❌")
else
    println("\n✅ Cell values remain reasonable")
end

# Save Julia result for further comparison
using MAT
matwrite("julia_step5_from_matlab_state.mat", Dict(
    "M" => M,
    "M_cell_713" => M[ih_track, jh_track, :],
    "t" => t_matlab + dt,
    "dt" => dt
))

println("\n✅ Saved Julia result to julia_step5_from_matlab_state.mat")

