"""
Simple test: Compare 3D (Nz=1) with archived 2D results

This test compares pre-computed results from 2D and 3D simulations
to validate backward compatibility without needing to run both simultaneously.
"""

using Test
using Printf
using JLD2

println("="^70)
println("Simple 2D vs 3D (Nz=1) Comparison")
println("="^70)

# Test parameters
const Nx = 20
const Ny = 20
const tmax = 0.01
const TEST_TOL_ABS = 1e-10
const TEST_TOL_REL = 1e-8

# Step 1: Generate 2D golden file (if needed)
golden_file_2d = joinpath(@__DIR__, "goldenfiles", "test_2d_nz1_comparison.jld2")

if !isfile(golden_file_2d)
    println("\n[folder] Generating 2D golden file...")
    println("   This requires running the archived 2D code separately")
    println("   Run from HyQMOM_2D_archive:")
    println("   julia --project=. -e 'using HyQMOM, JLD2; ...

'")
    @warn "Golden file not found: $golden_file_2d"
    @info "Skipping comparison - will just test that 3D with Nz=1 runs"
    HAS_GOLDEN = false
else
    println("\n[OK] Found 2D golden file")
    HAS_GOLDEN = true
end

# Step 2: Run 3D simulation with Nz=1
println("\n" * "="^70)
println("Running 3D simulation with Nz=1...")
println("="^70)

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using HyQMOM

results_3d = run_simulation(
    Nx=Nx,
    Ny=Ny,
    Nz=1,
    tmax=tmax,
    num_workers=1,
    verbose=false,
    homogeneous_z=true,
    enable_plots=false
)

println("\n3D (Nz=1) Results:")
println("  Final time: $(results_3d[:final_time])")
println("  Time steps: $(results_3d[:time_steps])")
println("  M shape: $(size(results_3d[:M]))")

@testset "3D with Nz=1 Basic Checks" begin
    @test results_3d[:M] !== nothing
    @test size(results_3d[:M], 1) == Np
    @test size(results_3d[:M], 2) == Np
    @test size(results_3d[:M], 3) == 1  # Nz=1
    @test size(results_3d[:M], 4) == 35
    @test results_3d[:final_time] ≈ tmax atol=1e-10
    @test all(results_3d[:M][:, :, 1, 1] .> 0)  # Density positive
    
    println("  [OK] All basic checks pass")
end

# Step 3: Compare if golden file exists
if HAS_GOLDEN
    println("\n" * "="^70)
    println("Comparing with 2D golden file...")
    println("="^70)
    
    # Load 2D golden data
    golden_data = load(golden_file_2d)
    M_2d = golden_data["M"]
    final_time_2d = golden_data["final_time"]
    time_steps_2d = golden_data["time_steps"]
    
    println("\n2D Golden Results:")
    println("  Final time: $final_time_2d")
    println("  Time steps: $time_steps_2d")
    println("  M shape: $(size(M_2d))")
    
    # Extract z-slice from 3D
    M_3d_slice = results_3d[:M][:, :, 1, :]
    
    @testset "2D vs 3D (Nz=1) Comparison" begin
        # Compare dimensions
        @test size(M_3d_slice) == size(M_2d)
        
        # Compare final times
        @test abs(results_3d[:final_time] - final_time_2d) < 1e-10
        
        # Compare moment arrays
        abs_diff = abs.(M_3d_slice .- M_2d)
        rel_diff = abs_diff ./ (abs.(M_2d) .+ 1e-30)
        
        max_abs_diff = maximum(abs_diff)
        max_rel_diff = maximum(rel_diff)
        
        idx_max = argmax(abs_diff)
        i_max, j_max, m_max = Tuple(idx_max)
        
        println("\nDifference Statistics:")
        println("  Max absolute difference: $max_abs_diff")
        println("  Max relative difference: $max_rel_diff")
        println("  Location: (i=$i_max, j=$j_max, moment=$m_max)")
        println("  2D value: $(M_2d[i_max, j_max, m_max])")
        println("  3D value: $(M_3d_slice[i_max, j_max, m_max])")
        
        @test max_abs_diff < TEST_TOL_ABS || max_rel_diff < TEST_TOL_REL
        
        if max_abs_diff < TEST_TOL_ABS && max_rel_diff < TEST_TOL_REL
            println("\n[OK] 3D (Nz=1) MATCHES 2D within tolerance!")
        else
            println("\n[WARNING]  Differences found")
        end
    end
else
    println("\n" * "="^70)
    println("Skipping comparison (no golden file)")
    println("To generate golden file, run:")
    println("  cd HyQMOM_2D_archive")
    println("  julia --project=. -e 'include(\"test/create_2d_golden.jl\")'")
    println("="^70)
end

println("\n" * "="^70)
println("Test Complete!")
println("="^70)

