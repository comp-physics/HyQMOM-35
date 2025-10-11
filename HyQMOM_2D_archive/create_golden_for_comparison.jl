"""
Generate golden file from 2D archived code for comparison with 3D (Nz=1)
"""

println("="^70)
println("Generating 2D Golden File for Comparison")
println("="^70)

using JLD2
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
using HyQMOM

# Test parameters - must match test_2d_vs_3d_simple.jl
const Np = 20
const tmax = 0.01

println("\nRunning 2D simulation...")
println("  Np = $Np")
println("  tmax = $tmax")

results = run_simulation(
    Np=Np,
    tmax=tmax,
    num_workers=1,
    verbose=false,
    enable_plots=false
)

println("\nResults:")
println("  Final time: $(results[:final_time])")
println("  Time steps: $(results[:time_steps])")
println("  M shape: $(size(results[:M]))")

# Save to golden file
output_dir = joinpath(@__DIR__, "..", "HyQMOM.jl", "test", "goldenfiles")
mkpath(output_dir)
output_file = joinpath(output_dir, "test_2d_nz1_comparison.jld2")

println("\nSaving to: $output_file")
jldsave(output_file;
        M=results[:M],
        final_time=results[:final_time],
        time_steps=results[:time_steps],
        xm=results[:xm],
        ym=results[:ym],
        Np=Np,
        tmax=tmax)

println("✓ Golden file saved successfully!")
println("="^70)

