#!/usr/bin/env julia
# Test the eigenvalues6_hyperbolic_3D path where the bug occurs

using RodneyHQMOM
using Printf

# Test moment input (from the debug output at step 4, cell 7,13)
M_in = [3.343246e-02, -1.438692e-02, 4.035618e-02, -4.728524e-02, 1.426107e-01,
        -1.960026e-02, -1.087988e-02, -3.175909e-02, -7.557277e-02, 3.408086e-02,
        1.894046e-02, 9.549695e-02, 5.758093e-02, 1.189652e-01, 1.931774e-01,
        0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e-02,
        0.000000e+00, 3.000000e-02, 0.000000e+00, 0.000000e+00, 9.000000e-02,
        0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
        0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 3.000000e-02]

println("="^70)
println("TESTING eigenvalues6_hyperbolic_3D")
println("="^70)

@printf("\nInput moments:\n")
@printf("  M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", M_in[1:5]...)
@printf("  M[3] (M200) = %.6e\n", M_in[3])
@printf("  M[5] (M400) = %.6e\n", M_in[5])

println("\nCalling eigenvalues6_hyperbolic_3D for axis=1...")
vpxmin, vpxmax, M_out = eigenvalues6_hyperbolic_3D(M_in, 1, 0, 0.0)

@printf("\nOutput moments:\n")
@printf("  M_out[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n", M_out[1:5]...)
@printf("  M_out[3] (M200) = %.6e\n", M_out[3])
@printf("  M_out[5] (M400) = %.6e\n", M_out[5])

@printf("\nEigenvalues: vpxmin=%.6e, vpxmax=%.6e\n", vpxmin, vpxmax)

if abs(M_out[3] - M_in[5]) < 1e-10
    println("\n❌ BUG CONFIRMED: M_out[3] equals input M_in[5]!")
    println("   This means M200 got replaced with M400!")
elseif abs(M_out[3]) > 1e6
    @printf("\n❌ BUG CONFIRMED: M_out[3] exploded to %.6e!\n", M_out[3])
elseif abs(M_out[3] - M_in[3]) > 1e-3
    @printf("\n⚠️  WARNING: M_out[3] changed significantly!\n")
    @printf("   Input:  %.6e\n", M_in[3])
    @printf("   Output: %.6e\n", M_out[3])
    @printf("   Ratio:  %.2fx\n", M_out[3] / M_in[3])
else
    println("\n✅ No obvious bug detected")
end

