#!/usr/bin/env julia
# Trace the first ~1e-9 difference in the calculation chain
# Step-by-step logging of ALL intermediate values during first time step

using RodneyHQMOM
using Printf

println("===== TRACING FIRST DIFFERENCE - JULIA =====\n")

# Setup (same as actual simulation)
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

# Initial conditions (exactly as in simulation)
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

# Background state
Mr_bg = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002)

# Crossing jets states
Uc = Ma / sqrt(2.0)
Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002)
Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002)

# Jet bounds
Csize = floor(Int, 0.1 * Np)
Mint = div(Np, 2) + 1
Maxt = div(Np, 2) + 1 + Csize
Minb = div(Np, 2) - Csize
Maxb = div(Np, 2)

# Fill grid
for ii in 1:nx
    gi = ii
    for jj in 1:ny
        gj = jj
        
        Mr = Mr_bg
        
        if gi >= Minb && gi <= Maxb && gj >= Minb && gj <= Maxb
            Mr = Mb
        end
        
        if gi >= Mint && gi <= Maxt && gj >= Mint && gj <= Maxt
            Mr = Mt
        end
        
        M[ii + halo, jj + halo, :] = Mr
    end
end

println("Initial conditions set")
println("Tracking cell (8,9) - one of the cells with differences\n")

# Track cell (8,9)
track_i = 8
track_j = 9
ih = track_i + halo
jh = track_j + halo

println("=== INITIAL STATE ===")
M_initial = M[ih, jh, :]
@printf("Cell (%d,%d) moments:\n", track_i, track_j)
@printf("  M000 = %.15e\n", M_initial[1])
@printf("  M100 = %.15e\n", M_initial[2])
@printf("  M200 = %.15e\n", M_initial[3])
@printf("  M210 = %.15e\n", M_initial[8])
@printf("  M130 = %.15e\n", M_initial[14])
println()

# Compute first time step dt
println("=== COMPUTING TIME STEP ===")

# Compute eigenvalues for all cells
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

@printf("vmax = %.15e\n", vmax)
@printf("dt = %.15e\n", dt)
println()

@printf("=== FLUX COMPUTATION FOR CELL (%d,%d) ===\n", track_i, track_j)

# Get moments at cell and neighbors for flux computation
M_cell = M[ih, jh, :]
M_left = M[ih-1, jh, :]
M_right = M[ih+1, jh, :]
M_down = M[ih, jh-1, :]
M_up = M[ih, jh+1, :]

println("Before flux computation:")
@printf("  Center M210 = %.15e\n", M_cell[8])
@printf("  Left M210   = %.15e\n", M_left[8])
@printf("  Right M210  = %.15e\n", M_right[8])
println()

# Compute fluxes (x-direction)
_, _, Fx_cell, _ = Flux_closure35_and_realizable_3D(M_cell, flag2D, Ma)
_, _, Fx_left, _ = Flux_closure35_and_realizable_3D(M_left, flag2D, Ma)
_, _, Fx_right, _ = Flux_closure35_and_realizable_3D(M_right, flag2D, Ma)

println("X-direction fluxes:")
@printf("  Fx_left[8]  = %.15e\n", Fx_left[8])
@printf("  Fx_cell[8]  = %.15e\n", Fx_cell[8])
@printf("  Fx_right[8] = %.15e\n", Fx_right[8])
println()

# HLL flux at left interface
Fhll_left, vpl_left, vpr_left = flux_HLL(M_left, M_cell, vpxmin[track_i-1,track_j],
                                          vpxmax[track_i-1,track_j], vpxmin[track_i,track_j],
                                          vpxmax[track_i,track_j], Fx_left, Fx_cell, 1)

# HLL flux at right interface
Fhll_right, vpl_right, vpr_right = flux_HLL(M_cell, M_right, vpxmin[track_i,track_j],
                                              vpxmax[track_i,track_j], vpxmin[track_i+1,track_j],
                                              vpxmax[track_i+1,track_j], Fx_cell, Fx_right, 1)

println("HLL fluxes:")
@printf("  Fhll_left[8]  = %.15e\n", Fhll_left[8])
@printf("  Fhll_right[8] = %.15e\n", Fhll_right[8])
println()

# Apply flux update
Mnpx = M_cell - dt/dx * (Fhll_right - Fhll_left)

println("After X-flux update:")
@printf("  Mnpx[8] = %.15e\n", Mnpx[8])
@printf("  Change = %.15e\n", Mnpx[8] - M_cell[8])
println()

# Y-direction (similar)
_, Fy_cell, _, _ = Flux_closure35_and_realizable_3D(M_cell, flag2D, Ma)
_, Fy_down, _, _ = Flux_closure35_and_realizable_3D(M_down, flag2D, Ma)
_, Fy_up, _, _ = Flux_closure35_and_realizable_3D(M_up, flag2D, Ma)

Fhll_down, _, _ = flux_HLL(M_down, M_cell, vpymin[track_i,track_j-1],
                            vpymax[track_i,track_j-1], vpymin[track_i,track_j],
                            vpymax[track_i,track_j], Fy_down, Fy_cell, 2)

Fhll_up, _, _ = flux_HLL(M_cell, M_up, vpymin[track_i,track_j],
                          vpymax[track_i,track_j], vpymin[track_i,track_j+1],
                          vpymax[track_i,track_j+1], Fy_cell, Fy_up, 2)

Mnpy = M_cell - dt/dy * (Fhll_up - Fhll_down)

println("After Y-flux update:")
@printf("  Mnpy[8] = %.15e\n", Mnpy[8])
@printf("  Change = %.15e\n", Mnpy[8] - M_cell[8])
println()

# Combine fluxes
Mnp_combined = Mnpx + Mnpy - M_cell

println("After combining X and Y fluxes:")
@printf("  Mnp_combined[8] = %.15e\n", Mnp_combined[8])
@printf("  Change from initial = %.15e\n", Mnp_combined[8] - M_cell[8])
println()

# Apply realizability
println("=== REALIZABILITY CORRECTION ===")
_, _, _, Mr = Flux_closure35_and_realizable_3D(Mnp_combined, flag2D, Ma)
println("After Flux_closure35_and_realizable_3D:")
@printf("  Mr[8] = %.15e\n", Mr[8])
@printf("  Change = %.15e\n", Mr[8] - Mnp_combined[8])
println()

_, _, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
println("After eigenvalues6_hyperbolic_3D (x):")
@printf("  Mr[8] = %.15e\n", Mr[8])
println()

_, _, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
println("After eigenvalues6_hyperbolic_3D (y):")
@printf("  Mr[8] = %.15e\n", Mr[8])
println()

_, _, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
println("After second Flux_closure35_and_realizable_3D:")
@printf("  Mr[8] = %.15e\n", Mr[8])
println()

# Apply collision
println("=== COLLISION OPERATOR ===")
MMC = collision35(Mr, dt, Kn)
println("After collision35:")
@printf("  MMC[8] = %.15e\n", MMC[8])
@printf("  Change = %.15e\n", MMC[8] - Mr[8])
println()

println("=== FINAL STATE ===")
@printf("Cell (%d,%d) after one time step:\n", track_i, track_j)
@printf("  Initial M210 = %.15e\n", M_cell[8])
@printf("  Final M210   = %.15e\n", MMC[8])
@printf("  Total change = %.15e\n", MMC[8] - M_cell[8])
println()

println("===== END JULIA TRACE =====")

# Also check M130 (moment 14)
println("\n=== COMPARISON FOR M130 (moment 14) ===")
@printf("  Initial M130 = %.15e\n", M_cell[14])
@printf("  Final M130   = %.15e\n", MMC[14])
@printf("  Total change = %.15e\n", MMC[14] - M_cell[14])
println()
