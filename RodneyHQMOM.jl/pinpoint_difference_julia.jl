#!/usr/bin/env julia
# Pinpoint exactly where M210 and M130 become nonzero in Julia
# Track through realizability corrections

using RodneyHQMOM
using Printf

println("===== PINPOINTING DIFFERENCE - JULIA =====\n")

# Same setup as simulation
Np = 20
rho = 0.7458726689777664  # From cell (8,9) at step 1
umean = -0.13731136169066
vmean = 0.0
wmean = 0.0
C200 = 0.6829490116088
C020 = 0.7458726689777664
C002 = 0.7458726689777664
C110 = 0.0
C101 = 0.0
C011 = 0.0
Ma = 0.0
flag2D = 0

println("Testing InitializeM4_35 with actual simulation values:")
@printf("  rho = %.15e\n", rho)
@printf("  umean = %.15e\n", umean)
@printf("  C200 = %.15e\n", C200)
@printf("  C110 = %.15e (should be 0)\n\n", C110)

# Call InitializeM4_35
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002)

println("Output from InitializeM4_35:")
@printf("  M210 = %.15e\n", M[8])
@printf("  M130 = %.15e\n\n", M[14])

# Now apply realizability corrections
println("=== APPLYING REALIZABILITY ===")

# First Flux_closure35
_, _, _, Mr = Flux_closure35_and_realizable_3D(M, flag2D, Ma)
println("After Flux_closure35_and_realizable_3D:")
@printf("  M210 = %.15e\n", Mr[8])
@printf("  M130 = %.15e\n\n", Mr[14])

# Eigenvalues X
_, _, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
println("After eigenvalues6_hyperbolic_3D (x-direction):")
@printf("  M210 = %.15e\n", Mr[8])
@printf("  M130 = %.15e\n\n", Mr[14])

# Eigenvalues Y
_, _, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
println("After eigenvalues6_hyperbolic_3D (y-direction):")
@printf("  M210 = %.15e\n", Mr[8])
@printf("  M130 = %.15e\n\n", Mr[14])

# Second Flux_closure35
_, _, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
println("After second Flux_closure35_and_realizable_3D:")
@printf("  M210 = %.15e\n", Mr[8])
@printf("  M130 = %.15e\n\n", Mr[14])

# Apply collision
dt = 0.02886751345948128
Kn = 1.0
Mr_final = collision35(Mr, dt, Kn)

println("After collision35:")
@printf("  M210 = %.15e\n", Mr_final[8])
@printf("  M130 = %.15e\n\n", Mr_final[14])

println("=== SUMMARY ===")
@printf("Initial (from InitializeM4_35):  M210 = %.15e, M130 = %.15e\n", M[8], M[14])
@printf("Final (after all corrections):   M210 = %.15e, M130 = %.15e\n", Mr_final[8], Mr_final[14])

println("\n===== END JULIA =====")

# Extra diagnostics
println("\n=== DETAILED DIAGNOSTICS ===")
println("All moments from InitializeM4_35:")
for i in 1:length(M)
    if abs(M[i]) < 1e-12 && M[i] != 0.0
        @printf("  M[%2d] = %.15e  (tiny nonzero)\n", i, M[i])
    end
end
