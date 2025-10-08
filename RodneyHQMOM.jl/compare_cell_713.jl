#!/usr/bin/env julia
# Compare cell (7,13) after step 4

using RodneyHQMOM
using MAT
using Printf

println("Running Julia to t=0.074...")
result = run_simulation(Np=20, tmax=0.074, verbose=false)

M = result[:M]
t = result[:final_time]
steps = result[:time_steps]

println("\n=== Julia Cell (7,13) after step $steps ===")
@printf("Time: %.15e\n", t)
println("M[7,13,1:10]:")
for k = 1:10
    @printf("  M[%d] = %.15e\n", k, M[7,13,k])
end

# Compare with MATLAB
if isfile("../matlab_cell713.mat")
    println("\n" * "="^60)
    println("COMPARISON WITH MATLAB")
    println("="^60)
    
    matlab_data = matread("../matlab_cell713.mat")
    M_matlab = matlab_data["result"]["M"]
    t_matlab = matlab_data["result"]["t"]
    steps_matlab = Int(matlab_data["result"]["steps"])
    
    @printf("\nSteps: MATLAB=%d, Julia=%d\n", steps_matlab, steps)
    @printf("Time: MATLAB=%.15e, Julia=%.15e\n", t_matlab, t)
    
    println("\nMoment comparison at cell (7,13):")
    for k = 1:10
        m_val = M_matlab[7,13,k]
        j_val = M[7,13,k]
        diff = abs(m_val - j_val)
        
        if abs(m_val) > 1e-30
            rel_err = diff / abs(m_val) * 100
        else
            rel_err = 0.0
        end
        
        status = rel_err < 1e-10 ? "✅" : (rel_err < 1e-6 ? "⚠️ " : "❌")
        @printf("%s M[%2d]: MATLAB=%.15e, Julia=%.15e, Diff=%.2e (%.2e%%)\n",
                status, k, m_val, j_val, diff, rel_err)
    end
    
    # Check if any moment is unreasonably large
    max_m_matlab = maximum(abs.(M_matlab[7,13,:]))
    max_m_julia = maximum(abs.(M[7,13,:]))
    
    println("\n" * "="^60)
    @printf("Max moment magnitude: MATLAB=%.2e, Julia=%.2e\n", max_m_matlab, max_m_julia)
    
    if max_m_julia > 1e6
        println("⚠️  WARNING: Julia moments are unreasonably large!")
    elseif max_m_matlab > 1e6
        println("⚠️  WARNING: MATLAB moments are unreasonably large!")
    else
        println("✅ Moment magnitudes are reasonable")
    end
end
