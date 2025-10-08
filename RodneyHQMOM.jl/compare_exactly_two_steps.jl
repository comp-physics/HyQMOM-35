#!/usr/bin/env julia
# Force exactly TWO time steps for comparison

using RodneyHQMOM
using MAT
using Printf

println("Running Julia simulation for exactly TWO steps...")
result = run_simulation(Np=20, tmax=0.04, verbose=true)

M_final = result[:M]

println("\n=== Julia Results ===")
println("Steps completed: ", result[:time_steps])
@printf("Final time: %.15e\n", result[:final_time])

println("\nFinal M[10,10,1:5]:")
for i = 1:5
    @printf("  M[10,10,%d] = %.15e\n", i, M_final[10,10,i])
end

# Save Julia results
matwrite("julia_exactly_two_steps.mat", Dict(
    "final_moments" => M_final,
    "final_time" => result[:final_time],
    "time_steps" => result[:time_steps]
))

# Compare with MATLAB
if isfile("../matlab_exactly_two_steps.mat")
    println("\n" * "="^60)
    println("COMPARISON")
    println("="^60)
    
    matlab_data = matread("../matlab_exactly_two_steps.mat")
    matlab_M = matlab_data["result"]["final_moments"]
    matlab_t = matlab_data["result"]["final_time"]
    matlab_steps = Int(matlab_data["result"]["time_steps"])
    
    julia_M = M_final
    julia_t = result[:final_time]
    julia_steps = result[:time_steps]
    
    @printf("\nSteps: MATLAB=%d, Julia=%d\n", matlab_steps, julia_steps)
    @printf("Final time: MATLAB=%.15e, Julia=%.15e\n", matlab_t, julia_t)
    
    println("\nPoint (10,10) comparison:")
    for i = 1:5
        m_val = matlab_M[10,10,i]
        j_val = julia_M[10,10,i]
        diff = abs(m_val - j_val)
        rel_err = diff / (abs(m_val) + 1e-30) * 100
        
        status = rel_err < 1e-10 ? "✅" : (rel_err < 1e-6 ? "⚠️ " : "❌")
        @printf("%s M[%d]: MATLAB=%.15e, Julia=%.15e, Diff=%.2e (%.2e%%)\n", 
                status, i, m_val, j_val, diff, rel_err)
    end
    
    # Check for NaN
    if any(isnan.(julia_M))
        println("\n❌ Julia has NaN values!")
    elseif any(isinf.(julia_M))
        println("\n❌ Julia has Inf values!")
    else
        max_err = maximum(abs.(matlab_M .- julia_M) ./ (abs.(matlab_M) .+ 1e-30)) * 100
        if max_err < 1e-10
            println("\n✅ PERFECT MATCH to machine precision!")
        elseif max_err < 1e-6
            @printf("\n✅ GOOD: Max error %.2e%%\n", max_err)
        else
            @printf("\n⚠️  Max error %.2e%%\n", max_err)
        end
    end
end
