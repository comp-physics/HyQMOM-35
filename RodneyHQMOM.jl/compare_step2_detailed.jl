#!/usr/bin/env julia
# Detailed comparison of step 2 to find algorithmic mismatch

using RodneyHQMOM
using MAT
using Printf

println("="^60)
println("DETAILED STEP 2 COMPARISON")
println("="^60)

# Run Julia step 1
println("\nRunning Julia step 1...")
result1 = run_simulation(Np=20, tmax=0.021, verbose=false)
M_julia_step1 = result1[:M]
t1_julia = result1[:final_time]
steps1_julia = result1[:time_steps]

println("\n=== After Step 1 (Julia) ===")
@printf("Time: %.15e, Steps: %d\n", t1_julia, steps1_julia)
@printf("Sample cell (10,10): M(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n",
        M_julia_step1[10,10,1], M_julia_step1[10,10,2], M_julia_step1[10,10,3],
        M_julia_step1[10,10,4], M_julia_step1[10,10,5])

# Run Julia step 2
println("\nRunning Julia step 2...")
result2 = run_simulation(Np=20, tmax=0.04, verbose=false)
M_julia_step2 = result2[:M]
t2_julia = result2[:final_time]
steps2_julia = result2[:time_steps]

println("\n=== After Step 2 (Julia) ===")
@printf("Time: %.15e, Steps: %d\n", t2_julia, steps2_julia)
@printf("Sample cell (10,10): M(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n",
        M_julia_step2[10,10,1], M_julia_step2[10,10,2], M_julia_step2[10,10,3],
        M_julia_step2[10,10,4], M_julia_step2[10,10,5])

# Load MATLAB
if !isfile("../matlab_step2_detailed.mat")
    println("\n❌ ERROR: Run save_step2_detailed.m first")
    exit(1)
end

println("\nLoading MATLAB results...")
matlab_data = matread("../matlab_step2_detailed.mat")
M_matlab_step1 = matlab_data["result"]["M_step1"]
M_matlab_step2 = matlab_data["result"]["M_step2"]
t1_matlab = matlab_data["result"]["t1"]
t2_matlab = matlab_data["result"]["t2"]

println("\n=== After Step 1 (MATLAB) ===")
@printf("Time: %.15e\n", t1_matlab)
@printf("Sample cell (10,10): M(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n",
        M_matlab_step1[10,10,1], M_matlab_step1[10,10,2], M_matlab_step1[10,10,3],
        M_matlab_step1[10,10,4], M_matlab_step1[10,10,5])

println("\n=== After Step 2 (MATLAB) ===")
@printf("Time: %.15e\n", t2_matlab)
@printf("Sample cell (10,10): M(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n",
        M_matlab_step2[10,10,1], M_matlab_step2[10,10,2], M_matlab_step2[10,10,3],
        M_matlab_step2[10,10,4], M_matlab_step2[10,10,5])

# Compare step 1
println("\n" * "="^60)
println("STEP 1 COMPARISON")
println("="^60)

for k in 1:10
    diff_step1 = abs.(M_matlab_step1[:,:,k] .- M_julia_step1[:,:,k])
    max_diff = maximum(diff_step1)
    max_idx = argmax(diff_step1)
    
    m_val = M_matlab_step1[max_idx[1], max_idx[2], k]
    j_val = M_julia_step1[max_idx[1], max_idx[2], k]
    rel = abs(m_val) > 1e-15 ? max_diff / abs(m_val) * 100 : 0.0
    
    status = rel < 1e-10 ? "✅" : (rel < 1e-6 ? "✅" : "❌")
    @printf("%s M[%2d]: max_diff=%.2e, rel=%.2e%% at (%d,%d)\n",
            status, k, max_diff, rel, max_idx[1], max_idx[2])
end

# Compare step 2
println("\n" * "="^60)
println("STEP 2 COMPARISON")
println("="^60)

for k in 1:10
    diff_step2 = abs.(M_matlab_step2[:,:,k] .- M_julia_step2[:,:,k])
    max_diff = maximum(diff_step2)
    max_idx = argmax(diff_step2)
    
    m_val = M_matlab_step2[max_idx[1], max_idx[2], k]
    j_val = M_julia_step2[max_idx[1], max_idx[2], k]
    rel = abs(m_val) > 1e-15 ? max_diff / abs(m_val) * 100 : 0.0
    
    status = rel < 1e-10 ? "✅" : (rel < 1e-6 ? "✅" : (rel < 0.1 ? "⚠️ " : "❌"))
    @printf("%s M[%2d]: max_diff=%.2e, rel=%.2e%% at (%d,%d)\n",
            status, k, max_diff, rel, max_idx[1], max_idx[2])
    
    if rel > 1e-6
        @printf("    MATLAB: %.15e, Julia: %.15e\n", m_val, j_val)
    end
end

# Find cells with largest differences
println("\n" * "="^60)
println("CELLS WITH LARGEST DIFFERENCES (Step 2)")
println("="^60)

# Sum differences across all moments for each cell
total_diff = zeros(size(M_matlab_step2, 1), size(M_matlab_step2, 2))
for k in 1:35
    total_diff .+= abs.(M_matlab_step2[:,:,k] .- M_julia_step2[:,:,k])
end

# Find top 5 cells
top_indices = sortperm(vec(total_diff), rev=true)[1:5]
for (rank, idx_flat) in enumerate(top_indices)
    idx = CartesianIndices(total_diff)[idx_flat]
    i, j = idx[1], idx[2]
    
    println("\n#$rank: Cell ($i,$j), Total diff = $(total_diff[i,j])")
    println("  First 5 moments:")
    for k in 1:5
        m_val = M_matlab_step2[i,j,k]
        j_val = M_julia_step2[i,j,k]
        diff = abs(m_val - j_val)
        rel = abs(m_val) > 1e-15 ? diff / abs(m_val) * 100 : 0.0
        @printf("    M[%d]: MATLAB=%.6e, Julia=%.6e, diff=%.2e (%.2e%%)\n",
                k, m_val, j_val, diff, rel)
    end
end

println("\n" * "="^60)
println("ANALYSIS")
println("="^60)
println("If step 1 is machine precision but step 2 has larger errors,")
println("there's an algorithmic difference in how step 2 is computed.")
println("Check:")
println("  1. Collision operator (BGK)")
println("  2. Time step calculation (dt)")
println("  3. Flux computation order")
println("  4. Realizability correction logic")
