#!/usr/bin/env julia
# Trace ALL intermediate values during first time step
# Find exactly where MATLAB and Julia first diverge

using RodneyHQMOM
using Printf
using MAT

println("="^60)
println("TRACING INTERMEDIATE VALUES - JULIA")
println("="^60)
println()

# Setup
Np = 20
Kn = 1.0
Ma = 0.0
flag2D = 0
CFL = 0.5
dx = 2.0 / Np
dy = 2.0 / Np
Nmom = 35

# Initialize
halo = 1
nx = Np
ny = Np
M = zeros(nx+2*halo, ny+2*halo, Nmom)

# Initial conditions
U0 = 0.0; V0 = 0.0; W0 = 0.0
T = 1.0
rhol = 1.0
rhor = 0.01
r110 = 0.0
r101 = 0.0
r011 = 0.0

C200 = T
C020 = T
C002 = T
C110 = r110 * sqrt(C200 * C020)
C101 = r101 * sqrt(C200 * C002)
C011 = r011 * sqrt(C020 * C002)

println("=== INITIAL CONDITION PARAMETERS ===")
@printf("C200 = %.15e\n", C200)
@printf("C110 = %.15e\n", C110)
@printf("C101 = %.15e\n", C101)
@printf("C020 = %.15e\n", C020)
@printf("C011 = %.15e\n", C011)
@printf("C002 = %.15e\n\n", C002)

Mr_bg = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002)
Uc = Ma / sqrt(2)
Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002)
Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002)

println("=== INITIAL MOMENTS ===")
@printf("Mr_bg[8]  = %.15e\n", Mr_bg[8])
@printf("Mr_bg[14] = %.15e\n", Mr_bg[14])
@printf("Mt[8]  = %.15e\n", Mt[8])
@printf("Mt[14] = %.15e\n", Mt[14])
@printf("Mb[8]  = %.15e\n", Mb[8])
@printf("Mb[14] = %.15e\n\n", Mb[14])

# Fill grid
Csize = floor(Int, 0.1 * Np)
Mint = div(Np, 2) + 1
Maxt = div(Np, 2) + 1 + Csize
Minb = div(Np, 2) - Csize
Maxb = div(Np, 2)

for ii in 1:nx
    for jj in 1:ny
        Mr = Mr_bg
        if ii >= Minb && ii <= Maxb && jj >= Minb && jj <= Maxb
            Mr = Mb
        end
        if ii >= Mint && ii <= Maxt && jj >= Mint && jj <= Maxt
            Mr = Mt
        end
        M[ii + halo, jj + halo, :] = Mr
    end
end

# Track cell (8,9)
track_i = 8
track_j = 9
ih = track_i + halo
jh = track_j + halo

println("=== CELL ($track_i,$track_j) INITIAL STATE ===")
M_init = M[ih, jh, :]
@printf("M[8]  = %.15e\n", M_init[8])
@printf("M[14] = %.15e\n\n", M_init[14])

# Compute dt
vpxmin = zeros(nx, ny)
vpxmax = zeros(nx, ny)
vpymin = zeros(nx, ny)
vpymax = zeros(nx, ny)

println("=== COMPUTING EIGENVALUES (first few cells) ===")
for i in 1:min(3, nx)
    for j in 1:min(3, ny)
        ih_loop = i + halo
        jh_loop = j + halo
        MOM = M[ih_loop, jh_loop, :]
        
        _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
        vpxmin[i,j], vpxmax[i,j], _ = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
        vpymin[i,j], vpymax[i,j], _ = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
        
        if i == 1 && j == 1
            println("Cell (1,1):")
            @printf("  vpxmin = %.15e, vpxmax = %.15e\n", vpxmin[i,j], vpxmax[i,j])
            @printf("  vpymin = %.15e, vpymax = %.15e\n", vpymin[i,j], vpymax[i,j])
        end
    end
end

# Complete eigenvalue computation for all cells
for i in 1:nx
    for j in 1:ny
        ih_loop = i + halo
        jh_loop = j + halo
        MOM = M[ih_loop, jh_loop, :]
        _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
        vpxmin[i,j], vpxmax[i,j], _ = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
        vpymin[i,j], vpymax[i,j], _ = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
    end
end

vmax = maximum([abs.(vpxmin[:])..., abs.(vpxmax[:])..., abs.(vpymin[:])..., abs.(vpymax[:])...])
dt = CFL * min(dx, dy) / vmax

@printf("\nvmax = %.15e\n", vmax)
@printf("dt   = %.15e\n\n", dt)

# Flux computation
println("=== FLUX COMPUTATION (cell $track_i,$track_j) ===")
Mx = zeros(nx+2*halo, ny+2*halo, Nmom)
My = zeros(nx+2*halo, ny+2*halo, Nmom)
Fy = zeros(nx+2*halo, ny+2*halo, Nmom)
Fx = zeros(nx+2*halo, ny+2*halo, Nmom)

for i in 1:nx
    for j in 1:ny
        ih_loop = i + halo
        jh_loop = j + halo
        MOM = M[ih_loop, jh_loop, :]
        _, Fy_temp, Fx_temp, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
        Fy[ih_loop, jh_loop, :] = Fy_temp
        Fx[ih_loop, jh_loop, :] = Fx_temp
        Mx[ih_loop, jh_loop, :] = Mr
        My[ih_loop, jh_loop, :] = Mr
        
        if i == track_i && j == track_j
            println("After Flux_closure35:")
            @printf("  Mr[8]  = %.15e\n", Mr[8])
            @printf("  Mr[14] = %.15e\n", Mr[14])
            @printf("  Fx[8]  = %.15e\n", Fx_temp[8])
            @printf("  Fy[8]  = %.15e\n\n", Fy_temp[8])
        end
    end
end

# Halo exchange
M[:,1,:] = M[:,ny+1,:]
M[:,ny+2*halo,:] = M[:,halo+1,:]
M[1,:,:] = M[nx+1,:,:]
M[nx+2*halo,:,:] = M[halo+1,:,:]

# Flux updates
println("=== FLUX UPDATES ===")
Mnpx = pas_HLL(Mx[halo+1:halo+nx, halo+1:halo+ny, :],
               Fx[halo+1:halo+nx, halo+1:halo+ny, :],
               dt, dx, vpxmin, vpxmax, true, true)

Mnpy = pas_HLL(My[halo+1:halo+nx, halo+1:halo+ny, :],
               Fy[halo+1:halo+nx, halo+1:halo+ny, :],
               dt, dy, vpymin, vpymax, true, true)

println("Cell ($track_i,$track_j) after flux updates:")
@printf("  Mnpx[8] = %.15e\n", Mnpx[track_i, track_j, 8])
@printf("  Mnpy[8] = %.15e\n\n", Mnpy[track_i, track_j, 8])

# Combine
Mnp_temp = zeros(nx+2*halo, ny+2*halo, Nmom)
for i in 1:nx
    for j in 1:ny
        Mnp_temp[i+halo, j+halo, :] = Mnpx[i,j,:] + Mnpy[i,j,:] - M[i+halo, j+halo, :]
    end
end

println("After combining fluxes:")
@printf("  Mnp_temp[%d,%d,8] = %.15e\n\n", ih, jh, Mnp_temp[ih, jh, 8])

# Realizability
println("=== REALIZABILITY CORRECTIONS ===")
Mnp = zeros(nx+2*halo, ny+2*halo, Nmom)

for i in 1:nx
    for j in 1:ny
        ih_loop = i + halo
        jh_loop = j + halo
        MOM = Mnp_temp[ih_loop, jh_loop, :]
        
        _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
        _, _, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
        _, _, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
        _, _, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
        
        Mnp[ih_loop, jh_loop, :] = Mr
        
        if i == track_i && j == track_j
            println("Cell ($i,$j) after realizability:")
            @printf("  Mr[8]  = %.15e\n", Mr[8])
            @printf("  Mr[14] = %.15e\n\n", Mr[14])
        end
    end
end

# Halo exchange
Mnp[:,1,:] = Mnp[:,ny+1,:]
Mnp[:,ny+2*halo,:] = Mnp[:,halo+1,:]
Mnp[1,:,:] = Mnp[nx+1,:,:]
Mnp[nx+2*halo,:,:] = Mnp[halo+1,:,:]

# Bulk assignment
M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]

# Collision
println("=== COLLISION OPERATOR ===")
M_before_coll = M[ih, jh, :]
println("Before collision:")
@printf("  M[8]  = %.15e\n", M_before_coll[8])
@printf("  M[14] = %.15e\n", M_before_coll[14])

# Extract covariance inputs to collision
rho_coll = M_before_coll[1]
umean_coll = M_before_coll[2] / rho_coll
vmean_coll = M_before_coll[6] / rho_coll
wmean_coll = M_before_coll[16] / rho_coll
C200_coll = M_before_coll[3]/rho_coll - umean_coll^2
C020_coll = M_before_coll[10]/rho_coll - vmean_coll^2
C002_coll = M_before_coll[20]/rho_coll - wmean_coll^2
Theta_coll = (C200_coll + C020_coll + C002_coll) / 3

println("\nCovariance inputs to collision:")
@printf("  rho   = %.15e\n", rho_coll)
@printf("  umean = %.15e\n", umean_coll)
@printf("  C200  = %.15e\n", C200_coll)
@printf("  C020  = %.15e\n", C020_coll)
@printf("  C002  = %.15e\n", C002_coll)
@printf("  Theta = %.15e\n\n", Theta_coll)

MMC = collision35(M_before_coll, dt, Kn)

println("After collision:")
@printf("  M[8]  = %.15e\n", MMC[8])
@printf("  M[14] = %.15e\n\n", MMC[14])

# Save all intermediate values
matwrite("julia_intermediate_trace.mat", Dict(
    "M_init" => M_init,
    "M_before_coll" => M_before_coll,
    "MMC" => MMC,
    "dt" => dt,
    "vmax" => vmax,
    "rho_coll" => rho_coll,
    "umean_coll" => umean_coll,
    "C200_coll" => C200_coll,
    "C020_coll" => C020_coll,
    "C002_coll" => C002_coll,
    "Theta_coll" => Theta_coll,
    "track_i" => track_i,
    "track_j" => track_j
))

println("Saved to julia_intermediate_trace.mat")
println("\n" * "="^60)
println("END JULIA TRACE")
println("="^60)
