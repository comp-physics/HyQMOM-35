#!/usr/bin/env julia
# Minimal test case: Collision operator on isotropic Maxwellian
# This should produce exactly zero for cross-moments M210 and M130
# but numerical differences between Julia and MATLAB occur

using RodneyHQMOM
using LinearAlgebra
using Printf

println("===== MINIMAL COLLISION TEST - JULIA =====\n")

# Test Case 1: Isotropic Maxwellian (all correlations = 0)
println("Test Case 1: Isotropic Maxwellian")
println("-"^50)

# Initial conditions (isotropic, at rest)
rho = 1.0
umean = 0.0
vmean = 0.0
wmean = 0.0
Theta = 1.0  # Temperature

# Create initial Maxwellian with isotropic covariance
M_initial = InitializeM4_35(rho, umean, vmean, wmean, Theta, 0.0, 0.0, Theta, 0.0, Theta)

println("Initial moments (should be isotropic Maxwellian):")
@printf("  M000 = %.15e\n", M_initial[1])
@printf("  M100 = %.15e\n", M_initial[2])
@printf("  M200 = %.15e\n", M_initial[3])
@printf("  M210 = %.15e  <-- Cross-moment (should be ~0)\n", M_initial[8])
@printf("  M130 = %.15e  <-- Cross-moment (should be ~0)\n", M_initial[14])
println()

# Apply collision operator (should preserve Maxwellian)
dt = 0.01
Kn = 1.0
M_after = collision35(M_initial, dt, Kn)

println("After collision operator:")
@printf("  M000 = %.15e\n", M_after[1])
@printf("  M100 = %.15e\n", M_after[2])
@printf("  M200 = %.15e\n", M_after[3])
@printf("  M210 = %.15e  <-- Cross-moment\n", M_after[8])
@printf("  M130 = %.15e  <-- Cross-moment\n", M_after[14])
println()

println("Changes:")
@printf("  ΔM210 = %.15e\n", M_after[8] - M_initial[8])
@printf("  ΔM130 = %.15e\n", M_after[14] - M_initial[14])
println()

# Test Case 2: Apply S4toC4_3D_r directly
println("\nTest Case 2: Direct S4toC4_3D_r call")
println("-"^50)

# Standardized moments for Gaussian
S300=0.0; S210=0.0; S201=0.0; S120=0.0; S111=0.0; S102=0.0
S030=0.0; S021=0.0; S012=0.0; S003=0.0
S400=3.0; S310=0.0; S301=0.0; S220=1.0; S211=0.0; S202=1.0
S130=0.0; S121=0.0; S112=0.0; S103=0.0; S040=3.0; S031=0.0
S022=1.0; S013=0.0; S004=3.0

# Covariance (isotropic)
C200 = Theta
C020 = Theta
C002 = Theta
C110 = 0.0
C101 = 0.0
C011 = 0.0

println("Input covariance matrix:")
@printf("  C200 = %.15e, C110 = %.15e, C101 = %.15e\n", C200, C110, C101)
@printf("  C020 = %.15e, C011 = %.15e\n", C020, C011)
@printf("  C002 = %.15e\n", C002)
println()

# Call S4toC4_3D_r
C4 = S4toC4_3D_r(C200, C110, C101, C020, C011, C002,
                 S300, S210, S201, S120, S111, S102, S030, S021, S012, S003,
                 S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004)

# Extract specific central moments
C000, C100, C200_out, C300, C400,
C010, C110_out, C210, C310,
C020_out, C120, C220,
C030, C130,
C040,
C001, C101_out, C201, C301,
C002_out, C102, C202,
C003, C103,
C004,
C011_out, C111, C211,
C021, C121,
C031,
C012, C112,
C013,
C022 = M4_to_vars(C4)

println("Output central moments:")
@printf("  C110 = %.15e\n", C110_out)
@printf("  C210 = %.15e  <-- Should be 0\n", C210)
@printf("  C120 = %.15e  <-- Should be 0\n", C120)
@printf("  C130 = %.15e  <-- Should be 0\n", C130)
println()

# Test Case 3: Matrix square root
println("\nTest Case 3: Matrix square root behavior")
println("-"^50)

C2 = [C200 C110 C101; C110 C020 C011; C101 C011 C002]
println("Covariance matrix C2:")
display(C2)
println()

A = sqrt(C2)
println("Matrix square root A = sqrt(C2):")
display(A)
println()

println("Check: A*A (should equal C2):")
display(A * A)
println()

println("Difference: A*A - C2:")
display(A * A - C2)
println()

@printf("Max absolute error: %.15e\n", maximum(abs.(A * A - C2)))
println()

# Test Case 4: Compare sqrt behavior with complex input
println("\nTest Case 4: sqrt behavior with complex matrices")
println("-"^50)

# For a diagonal identity matrix, sqrt should be trivial
# But let's test with element-wise conversion
C2_complex = ComplexF64.(C2)
A_complex_elementwise = sqrt.(C2_complex)
println("Element-wise sqrt.(ComplexF64.(C2)):")
display(A_complex_elementwise)
println()

println("This is different from matrix sqrt!")
println("Matrix sqrt treats the whole matrix, not elements.")
println()

# Test Case 5: Small perturbation
println("\nTest Case 5: Effect of tiny perturbation")
println("-"^50)

# Add a tiny cross-correlation
epsilon = 1e-9
M_perturbed = InitializeM4_35(rho, umean, vmean, wmean, Theta, epsilon, 0.0, Theta, 0.0, Theta)

@printf("With tiny perturbation (r110 = %.2e):\n", epsilon)
@printf("  M210 = %.15e\n", M_perturbed[8])
@printf("  M130 = %.15e\n", M_perturbed[14])
println()

# Apply collision
M_perturbed_after = collision35(M_perturbed, dt, Kn)

println("After collision:")
@printf("  M210 = %.15e\n", M_perturbed_after[8])
@printf("  M130 = %.15e\n", M_perturbed_after[14])
println()

println("===== END JULIA TEST =====")

# Test Case 6: Direct comparison at problematic values
println("\nTest Case 6: Problematic moment values from actual simulation")
println("-"^50)

# These are the actual values from cell (8,9) at step 1
# where MATLAB gives ~1e-9 but Julia gives 0
rho_actual = 0.7458726689777664
umean_actual = -0.13731136169066
vmean_actual = 0.0  # Should be small
wmean_actual = 0.0
C200_actual = 0.6829490116088
C020_actual = 0.7458726689777664
C002_actual = 0.7458726689777664

M_actual = InitializeM4_35(rho_actual, umean_actual, vmean_actual, wmean_actual,
                            C200_actual, 0.0, 0.0, C020_actual, 0.0, C002_actual)

println("Actual simulation values:")
@printf("  M210 = %.15e\n", M_actual[8])
@printf("  M130 = %.15e\n", M_actual[14])
@printf("  (MATLAB gives M210 ≈ 1.8e-9, M130 ≈ 5.2e-10)\n")
println()
