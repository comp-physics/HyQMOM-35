#!/usr/bin/env julia
# Compare one time step between MATLAB and Julia

using RodneyHQMOM
using MAT
using Printf

println("Running Julia simulation for ONE step...")
result = run_simulation(Np=20, tmax=0.02, verbose=false)

println("\n=== Julia One Step Results ===")
println("Steps completed: ", result[:time_steps])
@printf("Final time: %.10e\n", result[:final_time])
avg_dt = result[:final_time] / result[:time_steps]
@printf("Average dt: %.10e\n", avg_dt)

M_final = result[:M]

println("\nFinal M[1,1,1:5]:")
for i = 1:5
    @printf("  M[1,1,%d] = %.15e\n", i, M_final[1,1,i])
end

println("\nFinal M[10,10,1:5]:")
for i = 1:5
    @printf("  M[10,10,%d] = %.15e\n", i, M_final[10,10,i])
end

println("\nFinal M[20,20,1:5]:")
for i = 1:5
    @printf("  M[20,20,%d] = %.15e\n", i, M_final[20,20,i])
end

# Save Julia results
println("\nSaving Julia results to julia_one_step.mat...")
matwrite("julia_one_step.mat", Dict(
    "final_moments" => M_final,
    "final_time" => result[:final_time],
    "time_steps" => result[:time_steps]
))

# Load and compare with MATLAB if available
if isfile("../matlab_one_step.mat")
    println("\n=== Comparing with MATLAB ===")
    matlab_data = matread("../matlab_one_step.mat")
    
    # Extract MATLAB results
    matlab_M = matlab_data["result"]["final_moments"]
    matlab_t = matlab_data["result"]["final_time"]
    matlab_steps = Int(matlab_data["result"]["time_steps"])
    matlab_dt = matlab_data["result"]["dt_history"]
    
    julia_M = M_final
    julia_t = result[:final_time]
    julia_steps = result[:time_steps]
    julia_dt = julia_t / julia_steps
    
    println("\nTime step comparison:")
    @printf("  MATLAB steps: %d, Julia steps: %d, Match: %s\n", 
            matlab_steps, julia_steps, matlab_steps == julia_steps ? "✅" : "❌")
    @printf("  MATLAB avg dt: %.15e\n", matlab_dt)
    @printf("  Julia avg dt:  %.15e\n", julia_dt)
    if matlab_steps == julia_steps
        @printf("  Difference: %.15e (%.2e%%)\n", abs(matlab_dt - julia_dt), 
                abs(matlab_dt - julia_dt)/matlab_dt * 100)
    else
        println("  ⚠️  Different number of steps - cannot compare dt directly")
    end
    
    println("\nFinal time comparison:")
    @printf("  MATLAB: %.15e\n", matlab_t)
    @printf("  Julia:  %.15e\n", julia_t)
    @printf("  Difference: %.15e (%.2e%%)\n", abs(matlab_t - julia_t), 
            abs(matlab_t - julia_t)/abs(matlab_t) * 100)
    
    println("\nMoment comparison at key points:")
    
    # Compare at (1,1)
    println("\nPoint (1,1):")
    for i = 1:5
        m_val = matlab_M[1,1,i]
        j_val = julia_M[1,1,i]
        diff = abs(m_val - j_val)
        rel_err = diff / (abs(m_val) + 1e-30) * 100
        @printf("  M[%d]: MATLAB=%.15e, Julia=%.15e, Diff=%.2e (%.2e%%)\n", 
                i, m_val, j_val, diff, rel_err)
    end
    
    # Compare at (10,10)
    println("\nPoint (10,10):")
    for i = 1:5
        m_val = matlab_M[10,10,i]
        j_val = julia_M[10,10,i]
        diff = abs(m_val - j_val)
        rel_err = diff / (abs(m_val) + 1e-30) * 100
        @printf("  M[%d]: MATLAB=%.15e, Julia=%.15e, Diff=%.2e (%.2e%%)\n", 
                i, m_val, j_val, diff, rel_err)
    end
    
    # Compare at (20,20)
    println("\nPoint (20,20):")
    for i = 1:5
        m_val = matlab_M[20,20,i]
        j_val = julia_M[20,20,i]
        diff = abs(m_val - j_val)
        rel_err = diff / (abs(m_val) + 1e-30) * 100
        @printf("  M[%d]: MATLAB=%.15e, Julia=%.15e, Diff=%.2e (%.2e%%)\n", 
                i, m_val, j_val, diff, rel_err)
    end
    
    # Overall statistics
    max_abs_diff = maximum(abs.(matlab_M .- julia_M))
    max_rel_err = maximum(abs.(matlab_M .- julia_M) ./ (abs.(matlab_M) .+ 1e-30)) * 100
    
    println("\n=== Overall Statistics ===")
    @printf("Max absolute difference: %.2e\n", max_abs_diff)
    @printf("Max relative error: %.2e%%\n", max_rel_err)
    
    if max_rel_err < 1e-10
        println("✅ EXCELLENT: Results match to machine precision!")
    elseif max_rel_err < 1e-6
        println("✅ GOOD: Results match to 6 decimal places")
    elseif max_rel_err < 1e-3
        println("⚠️  WARNING: Results differ at 3rd decimal place")
    else
        println("❌ ERROR: Significant differences detected")
    end
else
    println("\n⚠️  MATLAB results not found. Run compare_one_step.m first.")
end
