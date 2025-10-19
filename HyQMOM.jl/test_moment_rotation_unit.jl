#!/usr/bin/env julia
"""
Unit test: Verify moment rotation formulas are correct.

Test that rotate_moments_z90 followed by rotate_moments_z90_inverse
returns the original moment vector (within machine precision).
"""

using HyQMOM
using Printf

println("="^70)
println("UNIT TEST: Moment Rotation Inverse")
println("="^70)

# Create a test moment vector with non-zero values
M_orig = collect(1.0:35.0)  # Simple test: M[i] = i

println("\nOriginal moment vector (first 10):")
for i in 1:10
    @printf("  M[%2d] = %.1f\n", i, M_orig[i])
end

# Apply forward rotation
M_rot = HyQMOM.rotate_moments_z90(M_orig)

println("\nRotated moment vector (first 10):")
for i in 1:10
    @printf("  M_rot[%2d] = %.1f\n", i, M_rot[i])
end

# Apply inverse rotation
M_back = HyQMOM.rotate_moments_z90_inverse(M_rot)

println("\nRotated-back moment vector (first 10):")
for i in 1:10
    @printf("  M_back[%2d] = %.1f\n", i, M_back[i])
end

# Check error
max_error = maximum(abs.(M_back .- M_orig))
rms_error = sqrt(sum((M_back .- M_orig).^2) / length(M_orig))

println("\n" * "="^70)
println("RESULTS")
println("="^70)
@printf("Max error: %.6e\n", max_error)
@printf("RMS error: %.6e\n", rms_error)

if max_error < 1e-14
    println("\n✅ PASS: Rotation inverse is correct to machine precision!")
elseif max_error < 1e-10
    println("\n✓ PASS: Rotation inverse is numerically correct.")
else
    println("\n❌ FAIL: Rotation formulas have errors!")
    println("\nDetailed errors (showing non-zero errors > 1e-12):")
    for i in 1:35
        err = abs(M_back[i] - M_orig[i])
        if err > 1e-12
            @printf("  M[%2d]: %.6e (orig: %.1f, back: %.1f)\n", i, err, M_orig[i], M_back[i])
        end
    end
end

println("="^70)

