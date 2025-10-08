#!/usr/bin/env julia
# Debug realizability at step 4, cell (7,13)

using RodneyHQMOM
using MAT
using Printf

println("="^60)
println("REALIZABILITY DEBUG AT STEP 4, CELL (7,13)")
println("="^60)

# Run Julia to just before step 4
println("\nRunning Julia to t=0.073...")
result = run_simulation(Np=20, tmax=0.073, verbose=false)

M = result[:M]
t = result[:final_time]
steps = result[:time_steps]

println("\n=== Julia at t=$t (step $steps) ===")
println("Cell (7,13) M[1:5]:")
MOM_julia = M[7,13,:]
for k = 1:5
    @printf("  M[%d] = %.15e\n", k, MOM_julia[k])
end

# Test realizability on this cell
println("\n=== Testing Flux_closure35_and_realizable_3D ===")
@printf("Input M[3] = %.15e\n", MOM_julia[3])

flag2D = 0
Ma = 0.0
Mx, My, Mz, Mr = Flux_closure35_and_realizable_3D(MOM_julia, flag2D, Ma)

@printf("After realizability M[3] = %.15e\n", Mr[3])
@printf("Change: %.15e (%.2fx)\n", Mr[3] - MOM_julia[3], Mr[3] / MOM_julia[3])

# Compare with MATLAB
if isfile("../matlab_realizability_test.mat")
    println("\n" * "="^60)
    println("COMPARISON WITH MATLAB")
    println("="^60)
    
    matlab_data = matread("../matlab_realizability_test.mat")
    MOM_matlab = matlab_data["result"]["MOM_before"]
    Mr_matlab = matlab_data["result"]["Mr_after"]
    Mx_matlab = matlab_data["result"]["Mx"]
    My_matlab = matlab_data["result"]["My"]
    
    println("\nInput moments comparison:")
    for k = 1:5
        diff = abs(MOM_matlab[k] - MOM_julia[k])
        rel = abs(MOM_matlab[k]) > 1e-30 ? diff / abs(MOM_matlab[k]) * 100 : 0.0
        status = rel < 1e-6 ? "✅" : "❌"
        @printf("%s M[%d]: MATLAB=%.6e, Julia=%.6e, Diff=%.2e (%.2e%%)\n",
                status, k, MOM_matlab[k], MOM_julia[k], diff, rel)
    end
    
    println("\nOutput moments after realizability:")
    for k = 1:5
        diff = abs(Mr_matlab[k] - Mr[k])
        rel = abs(Mr_matlab[k]) > 1e-30 ? diff / abs(Mr_matlab[k]) * 100 : 0.0
        status = rel < 1.0 ? "✅" : (rel < 100.0 ? "⚠️ " : "❌")
        @printf("%s M[%d]: MATLAB=%.6e, Julia=%.6e, Diff=%.2e (%.2e%%)\n",
                status, k, Mr_matlab[k], Mr[k], diff, rel)
    end
    
    println("\n" * "="^60)
    println("KEY FINDING")
    println("="^60)
    
    @printf("M[3] change in MATLAB: %.6e → %.6e (%.2fx)\n", 
            MOM_matlab[3], Mr_matlab[3], Mr_matlab[3] / MOM_matlab[3])
    @printf("M[3] change in Julia:  %.6e → %.6e (%.2fx)\n",
            MOM_julia[3], Mr[3], Mr[3] / MOM_julia[3])
    
    if abs(Mr_matlab[3] - Mr[3]) / abs(Mr_matlab[3]) > 10.0
        println("\n❌ CRITICAL: Realizability produces VERY different results!")
        println("   This is the source of the moment explosion.")
    elseif abs(Mr_matlab[3] - Mr[3]) / abs(Mr_matlab[3]) > 1.0
        println("\n⚠️  WARNING: Realizability differs by >1%")
    else
        println("\n✅ Realizability produces similar results")
    end
end
