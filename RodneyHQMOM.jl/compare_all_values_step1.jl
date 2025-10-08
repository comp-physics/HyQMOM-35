#!/usr/bin/env julia
# Compare ALL values after step 1

using RodneyHQMOM
using MAT
using Printf

println("="^60)
println("COMPREHENSIVE COMPARISON - STEP 1")
println("="^60)

# Run Julia
println("\nRunning Julia for step 1...")
result_julia = run_simulation(Np=20, tmax=0.021, verbose=false)

M_julia = result_julia[:M]
t_julia = result_julia[:final_time]
steps_julia = result_julia[:time_steps]

# Extract derived quantities
println("Extracting derived quantities from Julia...")
Ny, Nx, Nmom = size(M_julia)

rho_julia = zeros(Ny, Nx)
u_julia = zeros(Ny, Nx)
v_julia = zeros(Ny, Nx)
w_julia = zeros(Ny, Nx)
T_julia = zeros(Ny, Nx)
C200_julia = zeros(Ny, Nx)
C020_julia = zeros(Ny, Nx)
C002_julia = zeros(Ny, Nx)

for i in 1:Ny
    for j in 1:Nx
        M = M_julia[i,j,:]
        
        # Basic quantities
        rho_julia[i,j] = M[1]  # M000
        u_julia[i,j] = M[2] / M[1]  # M100/M000
        v_julia[i,j] = M[6] / M[1]  # M010/M000
        w_julia[i,j] = M[16] / M[1]  # M001/M000
        
        # Central moments (variances)
        C200_julia[i,j] = M[3] - M[2]^2/M[1]  # M200 - M100^2/M000
        C020_julia[i,j] = M[10] - M[6]^2/M[1]  # M020 - M010^2/M000
        C002_julia[i,j] = M[20] - M[16]^2/M[1]  # M002 - M001^2/M000
        
        # Temperature
        T_julia[i,j] = (C200_julia[i,j] + C020_julia[i,j] + C002_julia[i,j]) / 3
    end
end

# Load MATLAB results
if !isfile("../matlab_step1_all_values.mat")
    println("\n❌ ERROR: MATLAB results not found. Run save_all_values_step1.m first.")
    exit(1)
end

println("Loading MATLAB results...")
matlab_data = matread("../matlab_step1_all_values.mat")

rho_matlab = matlab_data["result"]["rho"]
u_matlab = matlab_data["result"]["u"]
v_matlab = matlab_data["result"]["v"]
w_matlab = matlab_data["result"]["w"]
T_matlab = matlab_data["result"]["T"]
C200_matlab = matlab_data["result"]["C200"]
C020_matlab = matlab_data["result"]["C020"]
C002_matlab = matlab_data["result"]["C002"]
M_matlab = matlab_data["result"]["M"]
t_matlab = matlab_data["result"]["t"]
steps_matlab = Int(matlab_data["result"]["steps"])

println("\n" * "="^60)
println("BASIC COMPARISON")
println("="^60)
@printf("Steps: MATLAB=%d, Julia=%d %s\n", steps_matlab, steps_julia,
        steps_matlab == steps_julia ? "✅" : "❌")
@printf("Time: MATLAB=%.15e, Julia=%.15e\n", t_matlab, t_julia)

# Function to compute statistics
function compare_field(name, matlab_field, julia_field)
    diff = abs.(matlab_field .- julia_field)
    max_diff = maximum(diff)
    max_idx = argmax(diff)
    
    # Relative error
    rel_err = diff ./ (abs.(matlab_field) .+ 1e-30) * 100
    max_rel = maximum(rel_err)
    max_rel_idx = argmax(rel_err)
    
    println("\n" * "-"^60)
    println("$name:")
    @printf("  Max absolute diff: %.6e at (%d,%d)\n", max_diff, max_idx[1], max_idx[2])
    @printf("    MATLAB: %.15e\n", matlab_field[max_idx])
    @printf("    Julia:  %.15e\n", julia_field[max_idx])
    @printf("  Max relative error: %.6e%% at (%d,%d)\n", max_rel, max_rel_idx[1], max_rel_idx[2])
    @printf("    MATLAB: %.15e\n", matlab_field[max_rel_idx])
    @printf("    Julia:  %.15e\n", julia_field[max_rel_idx])
    
    # Status
    if max_rel < 1e-10
        println("  ✅ PERFECT: Match to machine precision")
        return true
    elseif max_rel < 1e-6
        println("  ✅ EXCELLENT: Match to 6+ decimals")
        return true
    elseif max_rel < 0.1
        println("  ⚠️  GOOD: Match to 3+ decimals")
        return true
    else
        println("  ❌ ERROR: Significant difference!")
        return false
    end
end

println("\n" * "="^60)
println("FIELD-BY-FIELD COMPARISON")
println("="^60)

all_good = true
all_good &= compare_field("Density (rho)", rho_matlab, rho_julia)
all_good &= compare_field("Velocity u", u_matlab, u_julia)
all_good &= compare_field("Velocity v", v_matlab, v_julia)
all_good &= compare_field("Velocity w", w_matlab, w_julia)
all_good &= compare_field("Temperature T", T_matlab, T_julia)
all_good &= compare_field("Variance C200", C200_matlab, C200_julia)
all_good &= compare_field("Variance C020", C020_matlab, C020_julia)
all_good &= compare_field("Variance C002", C002_matlab, C002_julia)

# Compare all moments
println("\n" * "="^60)
println("ALL MOMENTS COMPARISON")
println("="^60)

for k in 1:min(15, Nmom)
    M_matlab_k = M_matlab[:,:,k]
    M_julia_k = M_julia[:,:,k]
    
    diff = abs.(M_matlab_k .- M_julia_k)
    max_diff = maximum(diff)
    max_idx = argmax(diff)
    
    rel_err = diff ./ (abs.(M_matlab_k) .+ 1e-30) * 100
    max_rel = maximum(rel_err)
    
    status = max_rel < 1e-10 ? "✅" : (max_rel < 1e-6 ? "✅" : (max_rel < 0.1 ? "⚠️ " : "❌"))
    
    @printf("%s M[%2d]: max_diff=%.2e, max_rel=%.2e%% at (%d,%d)\n",
            status, k, max_diff, max_rel, max_idx[1], max_idx[2])
    
    if max_rel > 0.1
        @printf("    MATLAB: %.6e, Julia: %.6e\n", 
                M_matlab_k[max_idx], M_julia_k[max_idx])
        all_good = false
    end
end

# Sample point comparison
println("\n" * "="^60)
println("SAMPLE POINT (10,10) DETAILED COMPARISON")
println("="^60)

test_points = [(10,10), (7,13), (1,1), (20,20)]

for (i,j) in test_points
    println("\nCell ($i,$j):")
    @printf("  rho:  MATLAB=%.15e, Julia=%.15e, diff=%.2e\n",
            rho_matlab[i,j], rho_julia[i,j], abs(rho_matlab[i,j] - rho_julia[i,j]))
    @printf("  u:    MATLAB=%.15e, Julia=%.15e, diff=%.2e\n",
            u_matlab[i,j], u_julia[i,j], abs(u_matlab[i,j] - u_julia[i,j]))
    @printf("  v:    MATLAB=%.15e, Julia=%.15e, diff=%.2e\n",
            v_matlab[i,j], v_julia[i,j], abs(v_matlab[i,j] - v_julia[i,j]))
    @printf("  T:    MATLAB=%.15e, Julia=%.15e, diff=%.2e\n",
            T_matlab[i,j], T_julia[i,j], abs(T_matlab[i,j] - T_julia[i,j]))
    
    println("  First 5 moments:")
    for k in 1:5
        m_val = M_matlab[i,j,k]
        j_val = M_julia[i,j,k]
        diff = abs(m_val - j_val)
        rel = abs(m_val) > 1e-30 ? diff / abs(m_val) * 100 : 0.0
        status = rel < 1e-10 ? "✅" : (rel < 1e-6 ? "✅" : "❌")
        @printf("    %s M[%d]: MATLAB=%.6e, Julia=%.6e, diff=%.2e (%.2e%%)\n",
                status, k, m_val, j_val, diff, rel)
    end
end

println("\n" * "="^60)
println("FINAL VERDICT")
println("="^60)

if all_good
    println("✅ SUCCESS: All fields match within acceptable tolerance!")
    println("   Step 1 is working correctly.")
else
    println("❌ FAILURE: Significant differences detected in step 1!")
    println("   Need to investigate the diverging fields.")
end
