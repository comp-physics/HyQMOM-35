#!/usr/bin/env julia
# Debug step-by-step to find where NaN first appears

using RodneyHQMOM
using Printf

# Modify to print detailed info at each step
println("="^60)
println("STEP-BY-STEP DEBUG")
println("="^60)

# We'll run multiple simulations with increasing tmax to isolate the problem
test_times = [0.02, 0.04, 0.057, 0.074, 0.09, 0.095, 0.10]

for tmax in test_times
    println("\n" * "-"^60)
    @printf("Testing tmax = %.3f\n", tmax)
    println("-"^60)
    
    result = run_simulation(Np=20, tmax=tmax, verbose=false)
    
    steps = result[:time_steps]
    final_t = result[:final_time]
    M = result[:M]
    
    has_nan = any(isnan.(M))
    has_inf = any(isinf.(M))
    
    @printf("  Steps: %d\n", steps)
    @printf("  Final time: %.15e\n", final_t)
    @printf("  M[10,10,1]: %.15e\n", M[10,10,1])
    @printf("  M[10,10,2]: %.15e\n", M[10,10,2])
    @printf("  M[10,10,3]: %.15e\n", M[10,10,3])
    
    if has_nan
        println("  ❌ NaN detected!")
        nan_count = sum(isnan.(M))
        @printf("     NaN count: %d / %d\n", nan_count, length(M))
        break
    elseif has_inf
        println("  ❌ Inf detected!")
        break
    else
        println("  ✅ No NaN/Inf")
    end
end

println("\n" * "="^60)
println("Now let's compare dt values between MATLAB and Julia")
println("="^60)

# Run a simulation with detailed output
println("\nRunning Julia with verbose output...")
result = run_simulation(Np=20, tmax=0.09, verbose=true)
