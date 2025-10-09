#!/usr/bin/env julia
# Simple investigation: Why does Julia produce exact zeros?

using RodneyHQMOM
using LinearAlgebra
using Printf

println("="^60)
println("SIMPLE INVESTIGATION - Julia vs MATLAB zeros")
println("="^60)
println()

# Use exact values from cell (8,9) at step 1
rho = 0.7458726689777664
umean = -0.13731136169066
vmean = 0.0
wmean = 0.0
C200 = 0.6829490116088
C020 = 0.7458726689777664
C002 = 0.7458726689777664
C110 = 0.0  # Exactly zero - no correlation
C101 = 0.0
C011 = 0.0

println("=== INPUT PARAMETERS ===")
@printf("  rho   = %.15e\n", rho)
@printf("  umean = %.15e\n", umean)
@printf("  vmean = %.15e\n", vmean)
@printf("  wmean = %.15e\n", wmean)
@printf("  C200  = %.15e\n", C200)
@printf("  C110  = %.15e (zero)\n", C110)
@printf("  C101  = %.15e (zero)\n", C101)
@printf("  C020  = %.15e\n", C020)
@printf("  C011  = %.15e (zero)\n", C011)
@printf("  C002  = %.15e\n\n", C002)

# Call InitializeM4_35
println("=== CALLING InitializeM4_35 ===")
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002)

@printf("M size: %s\n", size(M))
@printf("M[1]  (M000) = %.15e\n", M[1])
@printf("M[2]  (M100) = %.15e\n", M[2])
@printf("M[3]  (M200) = %.15e\n", M[3])
@printf("M[6]  (M010) = %.15e\n", M[6])
@printf("M[8]  (M210) = %.15e  ← KEY MOMENT\n", M[8])
@printf("M[10] (M020) = %.15e\n", M[10])
@printf("M[14] (M130) = %.15e  ← KEY MOMENT\n", M[14])
@printf("M[16] (M001) = %.15e\n", M[16])
@printf("M[20] (M002) = %.15e\n\n", M[20])

# Check if exactly zero
println("=== PRECISION CHECK ===")
@printf("M[8] == 0.0?  %s\n", M[8] == 0.0)
@printf("M[14] == 0.0? %s\n", M[14] == 0.0)
@printf("abs(M[8]) < eps()?  %s (eps = %.3e)\n", abs(M[8]) < eps(), eps())
@printf("abs(M[14]) < eps()? %s\n\n", abs(M[14]) < eps())

# Matrix square root check
println("=== MATRIX SQUARE ROOT ===")
C2 = [C200 C110 C101; C110 C020 C011; C101 C011 C002]
println("Covariance matrix C2:")
display(C2)
println()

A = sqrt(C2)
println("\nsqrt(C2) = A:")
display(A)
println()

println("\nA * A:")
display(A * A)
println()

diff = A * A - C2
@printf("\nMax |A*A - C2| = %.15e\n\n", maximum(abs.(diff)))

# Eigenvalue check
println("=== EIGENVALUE ANALYSIS ===")
eigs = eigvals(C2)
println("Eigenvalues:")
display(eigs)
println()
@printf("\nAll eigenvalues > 0? %s\n", all(eigs .> 0))
@printf("Min eigenvalue: %.15e\n\n", minimum(eigs))

# BLAS configuration
println("=== BLAS/LAPACK INFO ===")
blas_config = LinearAlgebra.BLAS.get_config()
println("BLAS vendor: $(blas_config)")
println()

# Summary
println("="^60)
println("SUMMARY")
println("="^60)
@printf("Julia InitializeM4_35 produces:\n")
@printf("  M[8]  = %.15e  (MATLAB: should check)\n", M[8])
@printf("  M[14] = %.15e  (MATLAB: should check)\n", M[14])
println()
if M[8] == 0.0 && M[14] == 0.0
    println("✓ Both moments are EXACTLY zero")
    println("  This is mathematically correct for isotropic covariance")
    println("  But numerical residuals in MATLAB might produce ~1e-9")
else
    println("✗ Moments are NOT exactly zero")
    @printf("  M[8]  magnitude: %.3e\n", abs(M[8]))
    @printf("  M[14] magnitude: %.3e\n", abs(M[14]))
end
println("="^60)
