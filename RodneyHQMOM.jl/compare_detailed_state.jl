#!/usr/bin/env julia
# Compare detailed state between MATLAB and Julia

using RodneyHQMOM
using MAT
using Printf

println("="^60)
println("DETAILED STATE COMPARISON")
println("="^60)

# Run Julia
println("\nRunning Julia simulation...")
result_julia = run_simulation(Np=20, tmax=0.074, verbose=false)

M_julia = result_julia[:M]
t_julia = result_julia[:final_time]
steps_julia = result_julia[:time_steps]

# Load MATLAB
if !isfile("../matlab_detailed_state.mat")
    println("ERROR: MATLAB state file not found. Run save_detailed_state.m first.")
    exit(1)
end

println("Loading MATLAB results...")
matlab_data = matread("../matlab_detailed_state.mat")
M_matlab = matlab_data["result"]["M"]
t_matlab = matlab_data["result"]["t"]
steps_matlab = Int(matlab_data["result"]["steps"])

println("\n" * "="^60)
println("BASIC COMPARISON")
println("="^60)
@printf("Steps: MATLAB=%d, Julia=%d %s\n", steps_matlab, steps_julia, 
        steps_matlab == steps_julia ? "✅" : "❌")
@printf("Time: MATLAB=%.15e, Julia=%.15e\n", t_matlab, t_julia)

# Compare multiple cells
test_cells = [(1,1), (7,13), (10,10), (15,15), (20,20)]

for (i, j) in test_cells
    println("\n" * "="^60)
    @printf("CELL (%d,%d) COMPARISON\n", i, j)
    println("="^60)
    
    max_err = 0.0
    max_err_idx = 0
    
    for k = 1:min(15, size(M_matlab, 3))
        m_val = M_matlab[i,j,k]
        j_val = M_julia[i,j,k]
        diff = abs(m_val - j_val)
        
        if abs(m_val) > 1e-30
            rel_err = diff / abs(m_val) * 100
        else
            rel_err = diff > 1e-30 ? Inf : 0.0
        end
        
        if rel_err > max_err
            max_err = rel_err
            max_err_idx = k
        end
        
        # Only print if there's significant error
        if rel_err > 1e-6 || k <= 5
            status = rel_err < 1e-10 ? "✅" : (rel_err < 1.0 ? "⚠️ " : "❌")
            @printf("%s M[%2d]: MATLAB=%+.6e, Julia=%+.6e, Diff=%.2e (%.2e%%)\n",
                    status, k, m_val, j_val, diff, rel_err)
        end
    end
    
    if max_err < 1e-10
        println("→ ✅ EXCELLENT: All moments match to machine precision")
    elseif max_err < 1.0
        @printf("→ ⚠️  GOOD: Max error %.2e%% at M[%d]\n", max_err, max_err_idx)
    else
        @printf("→ ❌ ERROR: Max error %.2e%% at M[%d]\n", max_err, max_err_idx)
    end
end

# Overall statistics
println("\n" * "="^60)
println("OVERALL STATISTICS")
println("="^60)

# Find cells with largest errors
max_abs_diff = 0.0
max_rel_err = 0.0
max_abs_loc = (0,0,0)
max_rel_loc = (0,0,0)

for i = 1:size(M_matlab, 1)
    for j = 1:size(M_matlab, 2)
        for k = 1:min(15, size(M_matlab, 3))
            m_val = M_matlab[i,j,k]
            j_val = M_julia[i,j,k]
            diff = abs(m_val - j_val)
            
            if diff > max_abs_diff
                max_abs_diff = diff
                max_abs_loc = (i,j,k)
            end
            
            if abs(m_val) > 1e-30
                rel_err = diff / abs(m_val) * 100
                if rel_err > max_rel_err
                    max_rel_err = rel_err
                    max_rel_loc = (i,j,k)
                end
            end
        end
    end
end

@printf("Max absolute difference: %.2e at cell (%d,%d) moment %d\n", 
        max_abs_diff, max_abs_loc[1], max_abs_loc[2], max_abs_loc[3])
@printf("  MATLAB: %.6e, Julia: %.6e\n", 
        M_matlab[max_abs_loc...], M_julia[max_abs_loc...])

@printf("\nMax relative error: %.2e%% at cell (%d,%d) moment %d\n",
        max_rel_err, max_rel_loc[1], max_rel_loc[2], max_rel_loc[3])
@printf("  MATLAB: %.6e, Julia: %.6e\n",
        M_matlab[max_rel_loc...], M_julia[max_rel_loc...])

# Count cells with large errors
large_err_count = 0
for i = 1:size(M_matlab, 1)
    for j = 1:size(M_matlab, 2)
        for k = 1:min(15, size(M_matlab, 3))
            m_val = M_matlab[i,j,k]
            j_val = M_julia[i,j,k]
            if abs(m_val) > 1e-30
                rel_err = abs(m_val - j_val) / abs(m_val) * 100
                if rel_err > 10.0  # More than 10% error
                    large_err_count += 1
                end
            end
        end
    end
end

total_values = size(M_matlab, 1) * size(M_matlab, 2) * min(15, size(M_matlab, 3))
@printf("\nCells with >10%% error: %d / %d (%.2f%%)\n", 
        large_err_count, total_values, large_err_count / total_values * 100)

println("\n" * "="^60)
println("VERDICT")
println("="^60)

if max_rel_err < 1e-10
    println("✅ PERFECT: Results match to machine precision!")
elseif max_rel_err < 1.0
    println("✅ GOOD: Results match to within 1%")
elseif max_rel_err < 100.0
    println("⚠️  WARNING: Results differ by up to $(round(max_rel_err, digits=2))%")
else
    println("❌ CRITICAL: Results differ by $(round(max_rel_err, sigdigits=3))%")
    println("   Major divergence detected!")
end
