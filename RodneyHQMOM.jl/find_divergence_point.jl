#!/usr/bin/env julia
# Find exact point where Julia diverges from MATLAB using existing data

using RodneyHQMOM
using MAT
using Printf

println("="^70)
println("FINDING DIVERGENCE POINT: Comparing Julia vs MATLAB Step 5")
println("="^70)

# Load MATLAB step 5 data
println("\n📥 Loading MATLAB step 5 results...")
matlab_data = matread("../matlab_step5_all_values.mat")
M_matlab = matlab_data["result"]["M"]
t_matlab = matlab_data["result"]["t"]
steps_matlab = matlab_data["result"]["steps"]

println("  ✓ MATLAB: t=$(t_matlab), steps=$(steps_matlab)")
println("  ✓ M shape: ", size(M_matlab))

# Run Julia to approximately same point
println("\n🚀 Running Julia simulation to t=$(t_matlab)...")
result_julia = run_simulation(Np=20, tmax=t_matlab, verbose=false, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5)

M_julia = result_julia[:M]
t_julia = result_julia[:final_time]
steps_julia = result_julia[:time_steps]

println("  ✓ Julia: t=$(t_julia), steps=$(steps_julia)")

# Check if Julia crashed
if M_julia === nothing
    println("\n❌ Julia returned nothing - simulation failed")
    exit(1)
end

if any(isnan, M_julia)
    println("\n❌ Julia has NaN values - checking where...")
    
    # Find first NaN
    for k in 1:size(M_julia, 3)
        if any(isnan, M_julia[:,:,k])
            idx = findfirst(isnan, M_julia[:,:,k])
            i, j = Tuple(idx)
            @printf("  First NaN at cell (%d,%d), moment %d\n", i, j, k)
            break
        end
    end
    
    println("\n⚠️  Julia crashed - cannot compare final states")
    println("Need to compare at earlier step before crash")
    exit(1)
end

println("\n✅ Julia completed without NaN")

# Compare all moment values
println("\n" * "="^70)
println("COMPARING ALL CELLS AND MOMENTS")
println("="^70)

Ny, Nx, Nmom = size(M_matlab)
max_diff = 0.0
max_diff_location = (0, 0, 0)
max_rel_diff = 0.0

diffs = zeros(Ny, Nx, Nmom)

for i in 1:Ny
    for j in 1:Nx
        for k in 1:Nmom
            diff = abs(M_matlab[i,j,k] - M_julia[i,j,k])
            diffs[i,j,k] = diff
            
            if diff > max_diff
                max_diff = diff
                max_diff_location = (i, j, k)
            end
            
            if abs(M_matlab[i,j,k]) > 1e-10
                rel_diff = diff / abs(M_matlab[i,j,k])
                if rel_diff > max_rel_diff
                    max_rel_diff = rel_diff
                end
            end
        end
    end
end

@printf("\n📊 Overall Statistics:\n")
@printf("  Max absolute difference: %.6e\n", max_diff)
@printf("  Location: cell (%d,%d), moment %d\n", max_diff_location...)
@printf("  Max relative difference: %.2f%%\n", max_rel_diff * 100)

if max_diff > 1e-6
    println("\n⚠️  SIGNIFICANT DIFFERENCES FOUND")
    
    # Find top 10 differences
    println("\n📍 Top 10 largest differences:")
    sorted_diffs = sort(vec(diffs), rev=true)
    
    count = 0
    for val in sorted_diffs
        if count >= 10
            break
        end
        
        idx = findfirst(x -> x == val, diffs)
        if idx !== nothing
            i, j, k = Tuple(idx)
            matlab_val = M_matlab[i,j,k]
            julia_val = M_julia[i,j,k]
            
            @printf("  %2d. Cell (%2d,%2d), M[%2d]: MATLAB=%.6e, Julia=%.6e, diff=%.6e\n",
                    count+1, i, j, k, matlab_val, julia_val, val)
            
            diffs[i,j,k] = 0  # Mark as used
            count += 1
        end
    end
    
    # Check cell (7,13) specifically
    println("\n📌 Cell (7,13) detailed comparison:")
    i, j = 7, 13
    @printf("  Moment | MATLAB          | Julia           | Difference\n")
    @printf("  -------|-----------------|-----------------|----------------\n")
    for k in 1:min(15, Nmom)
        matlab_val = M_matlab[i,j,k]
        julia_val = M_julia[i,j,k]
        diff = abs(matlab_val - julia_val)
        @printf("  M[%2d]  | %15.6e | %15.6e | %15.6e\n", k, matlab_val, julia_val, diff)
    end
    
elseif max_diff > 1e-12
    println("\n✅ Small differences (likely roundoff)")
    @printf("   Max difference: %.6e (acceptable)\n", max_diff)
else
    println("\n✅✅ PERFECT MATCH! (within machine precision)")
end

# Save comparison data
println("\n💾 Saving comparison data...")
matwrite("julia_matlab_comparison_step5.mat", Dict(
    "M_matlab" => M_matlab,
    "M_julia" => M_julia,
    "diffs" => diffs,
    "max_diff" => max_diff,
    "max_diff_location" => collect(max_diff_location),
    "t_matlab" => t_matlab,
    "t_julia" => t_julia
))

println("✅ Saved to julia_matlab_comparison_step5.mat")

