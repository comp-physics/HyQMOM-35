#!/usr/bin/env julia
# Simplified trace: What does the collision operator see?
# This traces through ONE cell to see where nonzero cross-moments appear

using RodneyHQMOM
using Printf

println("===== TRACE COLLISION INPUT - JULIA =====\n")

# Setup
Np = 20
Kn = 1.0
Ma = 0.0
flag2D = 0
CFL = 0.5
dx = 2.0 / Np
dy = 2.0 / Np
Nmom = 35

# Initialize grid with halo
halo = 1
nx = Np
ny = Np
M = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)

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

Mr_bg = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002)
Uc = Ma / sqrt(2.0)
Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002)
Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002)

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

track_i = 8
track_j = 9
ih = track_i + halo
jh = track_j + halo

@printf("Cell (%d,%d) initial state:\n", track_i, track_j)
@printf("  M210 = %.15e\n", M[ih, jh, 8])
@printf("  M130 = %.15e\n\n", M[ih, jh, 14])

# Compute dt
vpxmin = zeros(Float64, nx, ny)
vpxmax = zeros(Float64, nx, ny)
vpymin = zeros(Float64, nx, ny)
vpymax = zeros(Float64, nx, ny)

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

vmax = maximum([abs.(vpxmin); abs.(vpxmax); abs.(vpymin); abs.(vpymax)])
dt = CFL * min(dx, dy) / vmax
@printf("dt = %.15e\n\n", dt)

# Flux computation
Mx = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
My = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
Fy = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
Fx = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)

for i in 1:nx
    for j in 1:ny
        ih_loop = i + halo
        jh_loop = j + halo
        MOM = M[ih_loop, jh_loop, :]
        _, Fy[ih_loop,jh_loop,:], Fx[ih_loop,jh_loop,:], Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
        Mx[ih_loop, jh_loop, :] = Mr
        My[ih_loop, jh_loop, :] = Mr
    end
end

# Halo exchange (periodic)
M[:,1,:] = M[:,ny+1,:]
M[:,ny+2*halo,:] = M[:,halo+1,:]
M[1,:,:] = M[nx+1,:,:]
M[nx+2*halo,:,:] = M[halo+1,:,:]

# X-direction update
Mnpx = pas_HLL(Mx[halo+1:halo+nx, halo+1:halo+ny, :],
               Fx[halo+1:halo+nx, halo+1:halo+ny, :],
               dt, dx, vpxmin, vpxmax, true, true)

# Y-direction update
Mnpy = pas_HLL(My[halo+1:halo+nx, halo+1:halo+ny, :],
               Fy[halo+1:halo+nx, halo+1:halo+ny, :],
               dt, dy, vpymin, vpymax, true, true)

# Combine
Mnp_temp = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
for i in 1:nx
    for j in 1:ny
        Mnp_temp[i+halo, j+halo, :] = Mnpx[i,j,:] + Mnpy[i,j,:] - M[i+halo, j+halo, :]
    end
end

@printf("After flux update at cell (%d,%d):\n", track_i, track_j)
@printf("  M210 = %.15e\n", Mnp_temp[ih, jh, 8])
@printf("  M130 = %.15e\n\n", Mnp_temp[ih, jh, 14])

# Realizability
Mnp = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
v6xmin = zeros(Float64, nx, ny)
v6xmax = zeros(Float64, nx, ny)
v6ymin = zeros(Float64, nx, ny)
v6ymax = zeros(Float64, nx, ny)

for i in 1:nx
    for j in 1:ny
        ih_loop = i + halo
        jh_loop = j + halo
        MOM = Mnp_temp[ih_loop, jh_loop, :]
        
        _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
        v6xmin[i,j], v6xmax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
        v6ymin[i,j], v6ymax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
        _, _, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
        
        Mnp[ih_loop, jh_loop, :] = Mr
    end
end

@printf("After realizability at cell (%d,%d):\n", track_i, track_j)
@printf("  M210 = %.15e\n", Mnp[ih, jh, 8])
@printf("  M130 = %.15e\n\n", Mnp[ih, jh, 14])

# Halo exchange
Mnp[:,1,:] = Mnp[:,ny+1,:]
Mnp[:,ny+2*halo,:] = Mnp[:,halo+1,:]
Mnp[1,:,:] = Mnp[nx+1,:,:]
Mnp[nx+2*halo,:,:] = Mnp[halo+1,:,:]

# Collision
println("=== BEFORE COLLISION ===")
M_before_collision = Mnp[ih, jh, :]
@printf("  M210 = %.15e\n", M_before_collision[8])
@printf("  M130 = %.15e\n\n", M_before_collision[14])

MMC = collision35(M_before_collision, dt, Kn)

println("=== AFTER COLLISION ===")
@printf("  M210 = %.15e\n", MMC[8])
@printf("  M130 = %.15e\n\n", MMC[14])

println("=== FINAL ANSWER ===")
@printf("Cell (%d,%d) moments going into collision operator:\n", track_i, track_j)
@printf("  M210 = %.15e\n", M_before_collision[8])
@printf("  M130 = %.15e\n", M_before_collision[14])

println("\n===== END JULIA TRACE =====")

# Additional diagnostic: Check InitializeM4_35 directly
println("\n=== DIAGNOSTIC: Direct InitializeM4_35 test ===")
rho_test = M_before_collision[1]
u_test = M_before_collision[2] / rho_test
v_test = M_before_collision[6] / rho_test
w_test = M_before_collision[16] / rho_test

# Compute covariance from moments
C200_test = M_before_collision[3] / rho_test - u_test^2
C020_test = M_before_collision[10] / rho_test - v_test^2
C002_test = M_before_collision[20] / rho_test - w_test^2
C110_test = M_before_collision[7] / rho_test - u_test * v_test
C101_test = M_before_collision[17] / rho_test - u_test * w_test
C011_test = M_before_collision[26] / rho_test - v_test * w_test

@printf("Reconstructed parameters:\n")
@printf("  rho = %.15e\n", rho_test)
@printf("  u = %.15e\n", u_test)
@printf("  C110 = %.15e\n", C110_test)

# What does InitializeM4_35 give us?
M_reconstructed = InitializeM4_35(rho_test, u_test, v_test, w_test,
                                   C200_test, C110_test, C101_test,
                                   C020_test, C011_test, C002_test)
@printf("\nReconstructed M210 = %.15e\n", M_reconstructed[8])
@printf("Actual M210        = %.15e\n", M_before_collision[8])
@printf("Difference         = %.15e\n", M_reconstructed[8] - M_before_collision[8])
