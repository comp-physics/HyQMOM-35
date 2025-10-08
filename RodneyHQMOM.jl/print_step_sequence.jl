#!/usr/bin/env julia
# Print the exact step sequence for Julia

using RodneyHQMOM
using Printf

println("Julia step sequence for tmax=0.09:")
println("="^60)

result = run_simulation(Np=20, tmax=0.09, verbose=true)

println("\n" * "="^60)
@printf("Final: %d steps, t = %.15e\n", result[:time_steps], result[:final_time])

if any(isnan.(result[:M]))
    println("❌ Result contains NaN")
else
    println("✅ No NaN")
end
