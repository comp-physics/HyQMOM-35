#!/usr/bin/env julia
# Compare TWO time steps between MATLAB and Julia
# Find exact point of divergence

using RodneyHQMOM
using MAT
using Printf

println("Running Julia simulation for TWO steps...")
result = run_simulation(Np=20, tmax=0.1, verbose=true)

M_final = result[:M]

println("\n=== Julia Two Steps Results ===")
println("Steps completed: ", result[:time_steps])
@printf("Final time: %.15e\n", result[:final_time])

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
println("\nSaving Julia results to julia_two_steps.mat...")
matwrite("julia_two_steps.mat", Dict(
    "final_moments" => M_final,
    "final_time" => result[:final_time],
    "time_steps" => result[:time_steps]
))

# Load and compare with MATLAB if available
if isfile("../matlab_two_steps.mat")
    println("\n" * "="^60)
    println("COMPARING WITH MATLAB")
    println("="^60)
    
    matlab_data = matread("../matlab_two_steps.mat")
    
    # Extract MATLAB results
    matlab_M = matlab_data["result"]["final_moments"]
    matlab_t = matlab_data["result"]["final_time"]
    matlab_steps = Int(matlab_data["result"]["time_steps"])
    
    julia_M = M_final
    julia_t = result[:final_time]
    julia_steps = result[:time_steps]
    
    println("\n📊 STEP COUNT COMPARISON")
    @printf("  MATLAB steps: %d\n", matlab_steps)
    @printf("  Julia steps:  %d\n", julia_steps)
    if matlab_steps == julia_steps
        println("  ✅ Same number of steps")
    else
        println("  ❌ DIFFERENT number of steps!")
        println("  This is a critical difference - time stepping differs")
    end
    
    println("\n⏱️  FINAL TIME COMPARISON")
    @printf("  MATLAB: %.15e\n", matlab_t)
    @printf("  Julia:  %.15e\n", julia_t)
    @printf("  Difference: %.15e\n", abs(matlab_t - julia_t))
    
    # Detailed moment comparison at multiple points
    println("\n" * "="^60)
    println("DETAILED MOMENT COMPARISON")
    println("="^60)
    
    test_points = [(1,1), (5,5), (10,10), (15,15), (20,20)]
    
    for (i, j) in test_points
        println("\n📍 Point ($i,$j):")
        max_err = 0.0
        max_err_idx = 0
        
        for k = 1:5
            m_val = matlab_M[i,j,k]
            j_val = julia_M[i,j,k]
            diff = abs(m_val - j_val)
            
            # Handle division by zero
            if abs(m_val) > 1e-30
                rel_err = diff / abs(m_val) * 100
            else
                rel_err = 0.0
            end
            
            if rel_err > max_err
                max_err = rel_err
                max_err_idx = k
            end
            
            # Color code based on error
            if rel_err < 1e-10
                status = "✅"
            elseif rel_err < 1e-6
                status = "⚠️ "
            else
                status = "❌"
            end
            
            @printf("  %s M[%d]: MATLAB=%.15e, Julia=%.15e, Diff=%.2e (%.2e%%)\n", 
                    status, k, m_val, j_val, diff, rel_err)
        end
        
        if max_err < 1e-10
            println("  → EXCELLENT: All moments match to machine precision")
        elseif max_err < 1e-6
            @printf("  → GOOD: Max error %.2e%% at M[%d]\n", max_err, max_err_idx)
        else
            @printf("  → ⚠️  WARNING: Max error %.2e%% at M[%d]\n", max_err, max_err_idx)
        end
    end
    
    # Overall statistics
    println("\n" * "="^60)
    println("OVERALL STATISTICS")
    println("="^60)
    
    max_abs_diff = maximum(abs.(matlab_M .- julia_M))
    
    # Compute relative error carefully
    rel_errors = abs.(matlab_M .- julia_M) ./ (abs.(matlab_M) .+ 1e-30)
    max_rel_err = maximum(rel_errors) * 100
    
    # Find location of max error
    max_idx = argmax(abs.(matlab_M .- julia_M))
    
    @printf("Max absolute difference: %.2e\n", max_abs_diff)
    @printf("Max relative error: %.2e%%\n", max_rel_err)
    @printf("Location of max error: (%d,%d,%d)\n", max_idx[1], max_idx[2], max_idx[3])
    @printf("  MATLAB value: %.15e\n", matlab_M[max_idx])
    @printf("  Julia value:  %.15e\n", julia_M[max_idx])
    
    println("\n" * "="^60)
    println("VERDICT")
    println("="^60)
    
    if max_rel_err < 1e-10
        println("✅ EXCELLENT: Results match to machine precision!")
        println("   All functions are working identically.")
    elseif max_rel_err < 1e-6
        println("✅ GOOD: Results match to 6 decimal places")
        println("   Minor numerical differences, likely acceptable.")
    elseif max_rel_err < 1e-3
        println("⚠️  WARNING: Results differ at 3rd decimal place")
        println("   Investigate the point of maximum error.")
    else
        println("❌ ERROR: Significant differences detected")
        println("   Major divergence - need detailed step-by-step analysis.")
    end
    
    # Check for NaN/Inf
    if any(isnan.(julia_M)) || any(isinf.(julia_M))
        println("\n⚠️  Julia results contain NaN or Inf!")
        nan_count = sum(isnan.(julia_M))
        inf_count = sum(isinf.(julia_M))
        @printf("   NaN count: %d\n", nan_count)
        @printf("   Inf count: %d\n", inf_count)
    end
    
    if any(isnan.(matlab_M)) || any(isinf.(matlab_M))
        println("\n⚠️  MATLAB results contain NaN or Inf!")
        nan_count = sum(isnan.(matlab_M))
        inf_count = sum(isinf.(matlab_M))
        @printf("   NaN count: %d\n", nan_count)
        @printf("   Inf count: %d\n", inf_count)
    end
    
else
    println("\n⚠️  MATLAB results not found. Run compare_two_steps.m first.")
end
