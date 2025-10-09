#!/usr/bin/env julia
# Investigate why Julia produces exact zeros where MATLAB produces ~1e-9
# Focus on InitializeM4_35 and S4toC4_3D_r

using RodneyHQMOM
using LinearAlgebra
using Printf

println("="^60)
println("INVESTIGATING ZERO vs 1e-9 - JULIA")
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
println("Covariance matrix inputs:")
@printf("  C200 = %.15e\n", C200)
@printf("  C110 = %.15e (exactly zero)\n", C110)
@printf("  C101 = %.15e (exactly zero)\n", C101)
@printf("  C020 = %.15e\n", C020)
@printf("  C011 = %.15e (exactly zero)\n", C011)
@printf("  C002 = %.15e\n\n", C002)

# Step 1: Call S4toC4_3D_r directly
println("=== STEP 1: S4toC4_3D_r ===")

# Standardized moments for Gaussian
S300=0.0; S210=0.0; S201=0.0; S120=0.0; S111=0.0; S102=0.0
S030=0.0; S021=0.0; S012=0.0; S003=0.0
S400=3.0; S310=0.0; S301=0.0; S220=1.0; S211=0.0; S202=1.0
S130=0.0; S121=0.0; S112=0.0; S103=0.0; S040=3.0; S031=0.0
S022=1.0; S013=0.0; S004=3.0

println("Calling S4toC4_3D_r...")
C4_output = S4toC4_3D_r(C200, C110, C101, C020, C011, C002,
                        S300, S210, S201, S120, S111, S102, S030, S021, S012, S003,
                        S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004)

@printf("C4_output size: %s\n", size(C4_output))
@printf("C4_output is a 5x5x5 array\n")
# Index 8 in flat MATLAB indexing corresponds to...
# Need to convert to 3D indices
@printf("C4_output[2,2,1] (M210 in MATLAB linearization) = %.15e\n", C4_output[2,2,1])
println()

# Step 2: Look at matrix square root in S4toC4_3D_r
println("=== STEP 2: Matrix Square Root Analysis ===")

# Build the covariance matrix
C2 = [C200 C110 C101; C110 C020 C011; C101 C011 C002]
println("Covariance matrix C2:")
display(C2)
println()

println("\nComputing sqrt(C2)...")
A = sqrt(C2)
println("Matrix square root A:")
display(A)
println()

println("\nChecking A*A - C2:")
diff = A*A - C2
display(diff)
@printf("Max difference: %.15e\n\n", maximum(abs.(diff)))

# Step 3: Check eigenvalues and eigenvectors
println("=== STEP 3: Eigenvalue Analysis ===")
eigenvals = eigvals(C2)
println("Eigenvalues of C2:")
display(eigenvals)
println()

eigvecs = eigvecs(C2)
println("\nEigenvectors of C2:")
display(eigvecs)
println()

# Step 4: Trace through InitializeM4_35
println("\n=== STEP 4: InitializeM4_35 Output ===")
M_output = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002)

@printf("M_output size: %s\n", size(M_output))
@printf("M_output[8]  = %.15e (M210)\n", M_output[8])
@printf("M_output[14] = %.15e (M130)\n\n", M_output[14])

# Step 5: Check if it's a precision/tolerance issue
println("=== STEP 5: Precision Check ===")
@printf("eps() = %.15e\n", eps())
@printf("floatmin() = %.15e\n", floatmin())
@printf("Is M[8] exactly zero? %s\n", M_output[8] == 0.0)
@printf("Is M[8] < eps? %s\n", abs(M_output[8]) < eps())
@printf("Is M[8] < 1e-12? %s\n\n", abs(M_output[8]) < 1e-12)

# Step 6: Manual calculation to understand the computation
println("=== STEP 6: Manual Calculation ===")
println("From central to raw moments:")

println("Looking for C210 (central moment)...")

# The formula involves sqrt(C200) * sqrt(C020) terms
@printf("sqrt(C200) = %.15e\n", sqrt(C200))
@printf("sqrt(C020) = %.15e\n", sqrt(C020))
@printf("sqrt(C002) = %.15e\n", sqrt(C002))

# Extract from C4_output
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
C022 = M4_to_vars(C4_output)

println("\nExtracted central moments:")
@printf("  C210 = %.15e\n", C210)
@printf("  C120 = %.15e\n", C120)
@printf("  C130 = %.15e\n", C130)

println("\nComputing raw moment M210 from central C210:")
println("  M210 = C210 + umean*C020 + 2*vmean*C110 + vmean^2*C100")
println("       + umean^2*C010 + 2*umean*vmean*C010 + vmean^2*umean*M000")

M210_manual = C210 + umean*C020_out + 2*vmean*C110_out
println("  Simplified (vmean=0): M210 = C210 + umean*C020")
@printf("  M210_manual = %.15e + %.15e*%.15e\n", C210, umean, C020_out)
@printf("  M210_manual = %.15e\n\n", M210_manual)

# Step 7: Check for BLAS/LAPACK differences
println("=== STEP 7: BLAS/LAPACK Configuration ===")
println("BLAS config:")
display(LinearAlgebra.BLAS.get_config())
println()

println("=== SUMMARY ===")
@printf("Julia produces M210 = %.15e\n", M_output[8])
@printf("This is %s zero\n", M_output[8] == 0.0 ? "EXACTLY" : "NOT")
@printf("Magnitude: %.3e\n", abs(M_output[8]))

println("\n" * "="^60)
println("END JULIA INVESTIGATION")
println("="^60)
