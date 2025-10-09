#!/usr/bin/env julia

using RodneyHQMOM
using MAT
using Printf

println("===== JULIA: Two Time Steps =====\n")

params = SimulationParams(
    Np = 20,
    Kn = 1.0,
    Ma = 0.0,
    tmax = 0.05,
    flag2D = 0,
    CFL = 0.5
)

# Run simulation
elapsed = @elapsed result = run_simulation(params)

println("Completed in $(@sprintf("%.2f", elapsed)) seconds")
println("Time steps: $(result[:time_steps])")
println("Final time: $(@sprintf("%.6f", result[:final_time]))")

# Get moments
M_julia = result[:M]

println("\nSample moments at cell (8,9):")
@printf("  M[8,9,1] = %.15e\n", M_julia[8,9,1])
@printf("  M[8,9,3] = %.15e\n", M_julia[8,9,3])
@printf("  M[8,9,8] = %.15e\n", M_julia[8,9,8])
@printf("  M[8,9,14] = %.15e\n", M_julia[8,9,14])

# Load MATLAB results and compare
if isfile("../matlab_two_steps.mat")
    println("\n=== COMPARING WITH MATLAB ===")
    
    matlab_data = matread("../matlab_two_steps.mat")
    M_matlab = matlab_data["M_matlab"]
    
    # Compare
    max_diff = maximum(abs.(M_julia - M_matlab))
    rel_diff = maximum(abs.(M_julia - M_matlab) ./ (abs.(M_matlab) .+ 1e-16))
    
    @printf("Max absolute difference: %.3e\n", max_diff)
    @printf("Max relative difference: %.3e\n", rel_diff)
    
    if max_diff < 1e-12
        println("\n✅ PERFECT MATCH! (machine precision)")
    elseif max_diff < 1e-6
        println("\n✅ EXCELLENT MATCH! (< 1e-6)")
    elseif max_diff < 1e-3
        println("\n⚠️  Good match but small differences (< 1e-3)")
    else
        println("\n❌ SIGNIFICANT DIFFERENCES!")
    end
    
    # Find cells with largest differences
    println("\nCells with largest differences:")
    diffs = abs.(M_julia - M_matlab)
    for i in 1:5
        idx = argmax(diffs)
        @printf("  Cell [%d,%d,%d]: MATLAB=%.6e, Julia=%.6e, diff=%.3e\n",
                idx[1], idx[2], idx[3], M_matlab[idx...], M_julia[idx...], diffs[idx...])
        diffs[idx...] = 0  # Remove to find next largest
    end
else
    println("\n⚠️  MATLAB results not found. Run MATLAB first.")
end

println("\n✅ Julia results computed")
