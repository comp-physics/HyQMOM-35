#!/usr/bin/env julia
# Test step-by-step progression to find where explosion occurs

using RodneyHQMOM
using Printf

println("="^60)
println("STEP-BY-STEP EXPLOSION TEST")
println("="^60)

# Test different tmax values
test_times = [0.021, 0.04, 0.06, 0.074, 0.09]

for tmax in test_times
    println("\n" * "-"^60)
    @printf("Testing tmax = %.3f\n", tmax)
    println("-"^60)
    
    try
        result = run_simulation(Np=20, tmax=tmax, verbose=false)
        
        M = result[:M]
        steps = result[:time_steps]
        t_final = result[:final_time]
        
        # Check for NaN/Inf
        has_nan = any(isnan, M)
        has_inf = any(isinf, M)
        
        if has_nan || has_inf
            println("❌ EXPLOSION DETECTED!")
            @printf("   Steps completed: %d, Final time: %.6f\n", steps, t_final)
            
            # Find which cells exploded
            for k in 1:min(10, size(M, 3))
                if any(isnan, M[:,:,k]) || any(isinf, M[:,:,k])
                    idx = findfirst(x -> isnan(x) || isinf(x), M[:,:,k])
                    if idx !== nothing
                        @printf("   M[%d] exploded at cell (%d,%d) = %.2e\n", 
                                k, idx[1], idx[2], M[idx[1], idx[2], k])
                    end
                end
            end
        else
            println("✅ SUCCESS")
            @printf("   Steps completed: %d, Final time: %.6f\n", steps, t_final)
            
            # Check max values
            max_M = maximum(abs.(M))
            @printf("   Max |M| = %.2e\n", max_M)
            
            # Check specific cells
            @printf("   Cell (7,13): M[1]=%.2e, M[3]=%.2e\n", M[7,13,1], M[7,13,3])
            @printf("   Cell (14,8): M[1]=%.2e, M[3]=%.2e\n", M[14,8,1], M[14,8,3])
        end
        
    catch e
        println("❌ ERROR during simulation!")
        println("   ", e)
    end
end

println("\n" * "="^60)
println("SUMMARY")
println("="^60)
println("If explosion occurs, it happens between the last successful")
println("tmax and the first failed tmax.")
