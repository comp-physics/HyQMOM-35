#!/usr/bin/env julia
# Compare Julia simulation against MATLAB golden files

using RodneyHQMOM
using MAT
using MPI
using Printf

println("="^70)
println("GOLDEN FILE COMPARISON: Julia vs MATLAB")
println("="^70)

# Initialize MPI
MPI.Init()

# Test 1: Single rank, Np=20
println("\n" * "="^70)
println("TEST 1: Single Rank (1x1), Np=20, tmax=0.1")
println("="^70)

golden_file = "../goldenfiles/goldenfile_mpi_1ranks_Np20_tmax100.mat"
if !isfile(golden_file)
    error("Golden file not found: $golden_file")
end

# Load MATLAB golden file
println("\n📥 Loading MATLAB golden file...")
matlab_data = matread(golden_file)
golden_data = matlab_data["golden_data"]
M_matlab = golden_data["moments"]["M"]
t_matlab = golden_data["parameters"]["final_time"]
steps_matlab = golden_data["parameters"]["time_steps"]

println("  ✓ MATLAB: t=$(t_matlab), steps=$(steps_matlab), M shape = ", size(M_matlab))

# Run Julia simulation
println("\n🚀 Running Julia simulation...")
result_julia = run_simulation(
    Np=20, 
    tmax=t_matlab, 
    num_workers=1,
    Kn=1.0, 
    Ma=0.0, 
    flag2D=0, 
    CFL=0.5,
    verbose=false
)

M_julia = result_julia[:M]
t_julia = result_julia[:final_time]
steps_julia = result_julia[:time_steps]

println("  ✓ Julia: t=$(t_julia), steps=$(steps_julia), M shape = ", size(M_julia))

# Compare results
println("\n📊 Comparison Results:")
println("-"^70)

# Check dimensions
if size(M_matlab) != size(M_julia)
    println("  ❌ DIMENSION MISMATCH!")
    println("     MATLAB: ", size(M_matlab))
    println("     Julia:  ", size(M_julia))
else
    println("  ✓ Dimensions match: ", size(M_julia))
end

# Check time
if abs(t_matlab - t_julia) > 1e-10
    println("  ⚠️  Time mismatch: MATLAB=$(t_matlab), Julia=$(t_julia)")
else
    println("  ✓ Time matches: t=$(t_julia)")
end

# Compute differences
diff = M_julia .- M_matlab
abs_diff = abs.(diff)
rel_diff = abs_diff ./ (abs.(M_matlab) .+ 1e-30)

max_abs_diff = maximum(abs_diff)
max_rel_diff = maximum(rel_diff)
mean_abs_diff = sum(abs_diff) / length(abs_diff)
mean_rel_diff = sum(rel_diff) / length(rel_diff)

# Find location of maximum difference
max_diff_idx = argmax(abs_diff)
max_diff_loc = Tuple(max_diff_idx)

println("\n  Absolute Differences:")
println("    Max:  ", @sprintf("%.6e", max_abs_diff), " at ", max_diff_loc)
println("    Mean: ", @sprintf("%.6e", mean_abs_diff))

println("\n  Relative Differences:")
println("    Max:  ", @sprintf("%.6e", max_rel_diff))
println("    Mean: ", @sprintf("%.6e", mean_rel_diff))

# Check for NaN or Inf
nan_count = sum(isnan, M_julia)
inf_count = sum(isinf, M_julia)

if nan_count > 0 || inf_count > 0
    println("\n  ❌ NUMERICAL ISSUES:")
    println("     NaN count: ", nan_count)
    println("     Inf count: ", inf_count)
else
    println("\n  ✓ No NaN or Inf values")
end

# Overall pass/fail
TOLERANCE_ABS = 1e-8
TOLERANCE_REL = 1e-6

println("\n" * "="^70)
if max_abs_diff < TOLERANCE_ABS && max_rel_diff < TOLERANCE_REL
    println("✅ TEST PASSED! Julia matches MATLAB within tolerance.")
    println("   Absolute tolerance: ", @sprintf("%.6e", TOLERANCE_ABS))
    println("   Relative tolerance: ", @sprintf("%.6e", TOLERANCE_REL))
else
    println("⚠️  TEST WARNING: Differences exceed tolerance")
    println("   Max abs diff: ", @sprintf("%.6e", max_abs_diff), 
            " (tolerance: ", @sprintf("%.6e", TOLERANCE_ABS), ")")
    println("   Max rel diff: ", @sprintf("%.6e", max_rel_diff),
            " (tolerance: ", @sprintf("%.6e", TOLERANCE_REL), ")")
    
    # Show details at max diff location
    println("\n   At max diff location ", max_diff_loc, ":")
    println("     MATLAB: ", @sprintf("%.10e", M_matlab[max_diff_idx]))
    println("     Julia:  ", @sprintf("%.10e", M_julia[max_diff_idx]))
    println("     Diff:   ", @sprintf("%.10e", diff[max_diff_idx]))
end
println("="^70)

# Detailed statistics by moment
println("\n📈 Detailed Statistics by Moment Component:")
println("-"^70)
nx, ny, nmom = size(M_julia)
for k in 1:min(10, nmom)  # Show first 10 moments
    moment_abs_diff = abs_diff[:, :, k]
    moment_rel_diff = rel_diff[:, :, k]
    
    max_abs = maximum(moment_abs_diff)
    max_rel = maximum(moment_rel_diff)
    mean_abs = sum(moment_abs_diff) / length(moment_abs_diff)
    
    @printf("  Moment %2d: max_abs=%.6e, max_rel=%.6e, mean_abs=%.6e\n",
            k, max_abs, max_rel, mean_abs)
end

println("\n" * "="^70)
println("Golden file comparison complete!")
println("="^70)

