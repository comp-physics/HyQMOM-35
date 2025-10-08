#!/usr/bin/env julia
# Debug by modifying simulation_runner to print every step

using RodneyHQMOM
using Printf

# We need to temporarily modify the code to print every step
# For now, let's just run with a modified symmetry_check_interval

println("Running with detailed step output...")
println("="^60)

# Run simulation with symmetry_check_interval=1 to print every step
result = run_simulation(Np=20, tmax=0.09, verbose=true, symmetry_check_interval=1)

println("\n" * "="^60)
println("RESULT:")
@printf("Steps: %d\n", result[:time_steps])
@printf("Final time: %.15e\n", result[:final_time])

if any(isnan.(result[:M]))
    println("❌ NaN detected")
else
    println("✅ No NaN")
end
