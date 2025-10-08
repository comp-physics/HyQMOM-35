#!/usr/bin/env julia
# Compare ONLY the first step in extreme detail

using RodneyHQMOM
using MAT
using Printf

println("="^60)
println("FIRST STEP ONLY - EXTREME DETAIL")
println("="^60)

# The first step goes to t=0.020412
# Let's compare at exactly that time

println("\nRunning Julia to t=0.020412...")
result = run_simulation(Np=20, tmax=0.020412, verbose=false)
M_julia = result[:M]
t_julia = result[:final_time]
steps_julia = result[:time_steps]

println("\n=== Julia Results ===")
@printf("Time: %.15e, Steps: %d\n", t_julia, steps_julia)

# Load MATLAB step 1 results
if !isfile("../matlab_step2_detailed.mat")
    println("\n❌ ERROR: Run save_step2_detailed.m first")
    exit(1)
end

matlab_data = matread("../matlab_step2_detailed.mat")
M_matlab = matlab_data["result"]["M_step1"]
t_matlab = matlab_data["result"]["t1"]

println("\n=== MATLAB Results ===")
@printf("Time: %.15e\n", t_matlab)

# Note: MATLAB ran to t=0.021 (2 steps), we need to compare at same time
# Let's run Julia to exactly 0.021
println("\nRunning Julia to t=0.021...")
result = run_simulation(Np=20, tmax=0.021, verbose=false)
M_julia = result[:M]
t_julia = result[:final_time]
steps_julia = result[:time_steps]

println("\n=== Julia Results (t=0.021) ===")
@printf("Time: %.15e, Steps: %d\n", t_julia, steps_julia)

println("\n" * "="^60)
println("MOMENT-BY-MOMENT COMPARISON AT t=0.021")
println("="^60)

# Compare all 35 moments
for k in 1:35
    diff = abs.(M_matlab[:,:,k] .- M_julia[:,:,k])
    max_diff = maximum(diff)
    max_idx = argmax(diff)
    
    m_val = M_matlab[max_idx[1], max_idx[2], k]
    j_val = M_julia[max_idx[1], max_idx[2], k]
    
    # Calculate relative error carefully
    if abs(m_val) > 1e-15
        rel = max_diff / abs(m_val) * 100
    else
        rel = 0.0
    end
    
    status = if max_diff < 1e-14
        "✅ PERFECT"
    elseif rel < 1e-10
        "✅ EXCELLENT"
    elseif rel < 1e-6
        "✅ GOOD"
    elseif rel < 0.01
        "⚠️  SMALL"
    else
        "❌ LARGE"
    end
    
    @printf("%s M[%2d]: max_diff=%.2e, rel=%.2e%% at (%d,%d)\n",
            status, k, max_diff, rel, max_idx[1], max_idx[2])
    
    if max_diff > 1e-10 || rel > 1e-6
        @printf("         MATLAB=%.15e, Julia=%.15e\n", m_val, j_val)
    end
end

# Check specific problematic cells
println("\n" * "="^60)
println("PROBLEMATIC CELLS ANALYSIS")
println("="^60)

problem_cells = [(8,9), (11,8), (10,9)]

for (i,j) in problem_cells
    println("\nCell ($i,$j):")
    println("  Moments with largest differences:")
    
    diffs = Float64[]
    for k in 1:35
        push!(diffs, abs(M_matlab[i,j,k] - M_julia[i,j,k]))
    end
    
    # Find top 5 differences
    top_k = sortperm(diffs, rev=true)[1:5]
    for k in top_k
        m_val = M_matlab[i,j,k]
        j_val = M_julia[i,j,k]
        diff = diffs[k]
        rel = abs(m_val) > 1e-15 ? diff / abs(m_val) * 100 : 0.0
        
        @printf("    M[%2d]: diff=%.2e (%.2e%%), MATLAB=%.6e, Julia=%.6e\n",
                k, diff, rel, m_val, j_val)
    end
end

println("\n" * "="^60)
println("CONCLUSION")
println("="^60)
println("If errors exist even at the FIRST step, the problem is in:")
println("  1. Initial conditions")
println("  2. First flux computation")
println("  3. First realizability correction")
println("  4. First collision step")
println("NOT in accumulation over multiple steps!")
