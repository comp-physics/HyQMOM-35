#!/usr/bin/env julia
# Test Julia code for 13 time steps

using RodneyHQMOM
using Printf

println("="^60)
println("Testing Julia Code for 13 Time Steps")
println("="^60)
println()

# Set up parameters
Np = 20
Kn = 1.0
Ma = 0.0
tmax = 0.12  # Should give ~12-13 steps
flag2D = 0
CFL = 0.5

println("Parameters:")
println("  Np = $Np")
println("  Kn = $Kn")
println("  Ma = $Ma")
println("  tmax = $tmax")
println("  CFL = $CFL")
println()

println("Running simulation...")
println()

elapsed = @elapsed begin
    result = run_simulation(
        Np = Np,
        Kn = Kn,
        Ma = Ma,
        tmax = tmax,
        flag2D = flag2D,
        CFL = CFL,
        num_workers = 1,
        verbose = true,
        save_output = false
    )
end

println()
println("="^60)
println("RESULTS")
println("="^60)
println()

println("✅ Simulation completed successfully!")
println()
println("Time:")
println("  Wall time: $(@sprintf("%.2f", elapsed)) seconds")
println("  Final time: $(@sprintf("%.6f", result[:final_time]))")
println("  Time steps: $(result[:time_steps])")
println()

# Check for NaN or Inf
M = result[:M]
has_nan = any(isnan, M)
has_inf = any(isinf, M)
max_moment = maximum(abs, M)

println("Moment statistics:")
println("  Has NaN: $(has_nan ? "❌ YES" : "✅ NO")")
println("  Has Inf: $(has_inf ? "❌ YES" : "✅ NO")")
@printf("  Max |M|: %.6e\n", max_moment)
println()

if max_moment > 1e6
    println("⚠️  WARNING: Large moment values detected")
    println()
elseif has_nan || has_inf
    println("❌ FAILED: NaN or Inf detected")
    println()
else
    println("✅ SUCCESS: All moments are finite and bounded!")
    println()
end

# Sample moments at a few cells
println("Sample moments at cell (8,9):")
@printf("  M[8,9,1]  = %.6e (density)\n", M[8,9,1])
@printf("  M[8,9,3]  = %.6e\n", M[8,9,3])
@printf("  M[8,9,8]  = %.6e (M210)\n", M[8,9,8])
@printf("  M[8,9,14] = %.6e (M130)\n", M[8,9,14])
println()

println("Sample moments at cell (10,10):")
@printf("  M[10,10,1]  = %.6e (density)\n", M[10,10,1])
@printf("  M[10,10,3]  = %.6e\n", M[10,10,3])
@printf("  M[10,10,8]  = %.6e (M210)\n", M[10,10,8])
@printf("  M[10,10,14] = %.6e (M130)\n", M[10,10,14])
println()

println("="^60)
println("Julia simulation completed $(result[:time_steps]) time steps")
println("="^60)
