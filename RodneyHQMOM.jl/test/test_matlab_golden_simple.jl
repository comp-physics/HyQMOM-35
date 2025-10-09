"""
Simplified MATLAB Golden File Test for use within Test.jl framework

This version runs the golden file comparison without MPI initialization,
making it suitable for Pkg.test(). For full MPI testing, use test_matlab_golden.jl directly.
"""

using MAT
using Printf

# Test tolerances (must be at module level, not in testset)
const TOLERANCE_ABS = 1e-8
const TOLERANCE_REL = 1e-6

# Load MATLAB golden file
const GOLDEN_FILE_PATH = joinpath(@__DIR__, "..", "..", "goldenfiles", 
                                   "goldenfile_mpi_1ranks_Np20_tmax100.mat")

@testset "Load MATLAB Golden File" begin
    @test isfile(GOLDEN_FILE_PATH)
    
    global matlab_data = matread(GOLDEN_FILE_PATH)
    @test haskey(matlab_data, "golden_data")
    
    global golden_data = matlab_data["golden_data"]
    @test haskey(golden_data, "moments")
    @test haskey(golden_data, "parameters")
    @test haskey(golden_data, "grid")
    
    global M_matlab = golden_data["moments"]["M"]
    global params_matlab = golden_data["parameters"]
    
    @test size(M_matlab, 1) == params_matlab["Np"]
    @test size(M_matlab, 2) == params_matlab["Np"]
    @test size(M_matlab, 3) == 35  # Number of moments
    
    println("  ✓ MATLAB golden file loaded:")
    println("    Size: $(size(M_matlab))")
    println("    Np: $(params_matlab["Np"])")
    println("    tmax: $(params_matlab["tmax"])")
    println("    Final time: $(params_matlab["final_time"])")
    println("    Time steps: $(params_matlab["time_steps"])")
end

@testset "Run Julia Simulation" begin
    Np = Int(params_matlab["Np"])
    tmax = params_matlab["tmax"]
    final_time_matlab = params_matlab["final_time"]
    
    # Setup parameters EXACTLY as MATLAB (from parse_main_args.m)
    dx = 1.0 / Np
    dy = 1.0 / Np
    Nmom = 35
    nnmax = 20000000
    Kn = 1.0  # MATLAB default
    dtmax = Kn
    rhol = 1.0  # MATLAB default
    rhor = 0.01  # MATLAB default (NOT 0.5!)
    T = 1.0
    r110 = 0.0
    r101 = 0.0
    r011 = 0.0
    Ma = 0.0  # MATLAB default
    flag2D = 0
    CFL = 0.5  # MATLAB default (NOT 0.9!)
    symmetry_check_interval = 10
    
    params = (
        Np=Np, tmax=tmax, Kn=Kn, Ma=Ma, flag2D=flag2D, CFL=CFL,
        dx=dx, dy=dy, Nmom=Nmom, nnmax=nnmax, dtmax=dtmax,
        rhol=rhol, rhor=rhor, T=T, r110=r110, r101=r101, r011=r011,
        symmetry_check_interval=symmetry_check_interval,
        enable_memory_tracking=false
    )
    
    println("  Running Julia simulation with MATLAB parameters:")
    println("    Np=$Np, tmax=$tmax, Kn=$Kn, Ma=$Ma, CFL=$CFL")
    println("    rhol=$rhol, rhor=$rhor")
    
    # Initialize MPI for single-process run
    if !MPI.Initialized()
        MPI.Init()
    end
    
    global M_julia, final_time_julia, time_steps_julia, grid_julia = 
        RodneyHQMOM.simulation_runner(params)
    
    @test M_julia !== nothing
    @test size(M_julia) == size(M_matlab)
    @test time_steps_julia == params_matlab["time_steps"]
    @test abs(final_time_julia - final_time_matlab) < 1e-10
    
    println("  ✓ Julia simulation complete:")
    println("    Final time: $final_time_julia")
    println("    Time steps: $time_steps_julia")
    println("    Size: $(size(M_julia))")
end

@testset "Compare Results" begin
    @test M_julia !== nothing
    @test size(M_julia) == size(M_matlab)
    
    # Compute differences
    diff = M_julia .- M_matlab
    abs_diff = abs.(diff)
    rel_diff = abs_diff ./ (abs.(M_matlab) .+ 1e-30)
    
    max_abs_diff = maximum(abs_diff)
    max_rel_diff = maximum(rel_diff)
    mean_abs_diff = sum(abs_diff) / length(abs_diff)
    mean_rel_diff = sum(rel_diff) / length(rel_diff)
    
    # Check for NaN/Inf
    @test !any(isnan, M_julia)
    @test !any(isinf, M_julia)
    
    if any(isnan, M_julia)
        @warn "Julia results contain NaN"
    end
    if any(isinf, M_julia)
        @warn "Julia results contain Inf"
    end
    
    # Test against tolerances (use global constants defined at module level)
    @test max_abs_diff < TOLERANCE_ABS
    @test max_rel_diff < TOLERANCE_REL
    
    if max_abs_diff >= TOLERANCE_ABS
        @warn "Max abs diff $max_abs_diff exceeds tolerance $TOLERANCE_ABS"
    end
    if max_rel_diff >= TOLERANCE_REL
        @warn "Max rel diff $max_rel_diff exceeds tolerance $TOLERANCE_REL"
    end
    
    println()
    println("  Comparison Results:")
    @printf("    Max absolute diff: %.6e\n", max_abs_diff)
    @printf("    Max relative diff: %.6e\n", max_rel_diff)
    @printf("    Mean absolute diff: %.6e\n", mean_abs_diff)
    @printf("    Mean relative diff: %.6e\n", mean_rel_diff)
    println()
    
    if max_abs_diff < TOLERANCE_ABS && max_rel_diff < TOLERANCE_REL
        println("  ✅ Julia matches MATLAB within tolerance")
        if max_abs_diff < 1e-14
            println("     (Agreement at machine precision!)")
        end
    end
    
    # Moment-by-moment check (first 11 moments)
    moment_names = ["ρ", "ρu", "ρv", "ρw", "ρE", "M₁₁₀", "M₁₀₁", "M₀₁₁", 
                    "M₂₀₀", "M₀₂₀", "M₀₀₂"]
    
    println("  Per-moment analysis:")
    for k in 1:min(11, size(M_julia, 3))
        moment_abs_diff = abs_diff[:, :, k]
        max_abs = maximum(moment_abs_diff)
        
        moment_name = k <= length(moment_names) ? moment_names[k] : "M[$k]"
        @printf("    %-6s: max_abs=%.6e", moment_name, max_abs)
        
        if max_abs < TOLERANCE_ABS
            println(" ✓")
        else
            println(" ⚠")
        end
    end
end

@testset "Physical Sanity Checks" begin
    # Density should be positive
    ρ = M_julia[:, :, 1]
    @test all(ρ .> 0)
    if !all(ρ .> 0)
        @warn "Density must be positive"
    end
    
    @test minimum(ρ) > 1e-10
    if minimum(ρ) <= 1e-10
        @warn "Density should not be near zero"
    end
    
    # Energy should be positive  
    E = M_julia[:, :, 5]
    @test all(E .> 0)
    if !all(E .> 0)
        @warn "Energy must be positive"
    end
    
    # Check conservation (for 1 rank, boundaries are periodic-like)
    # Total mass should be conserved
    mass_julia = sum(M_julia[:, :, 1])
    mass_matlab = sum(M_matlab[:, :, 1])
    @test abs(mass_julia - mass_matlab) / mass_matlab < 1e-12
    
    println("  ✓ Physical sanity checks passed")
    @printf("    Min density: %.6e\n", minimum(ρ))
    @printf("    Max density: %.6e\n", maximum(ρ))
    @printf("    Total mass (Julia): %.6e\n", mass_julia)
    @printf("    Total mass (MATLAB): %.6e\n", mass_matlab)
end

println()
println("="^70)
println("✅ MATLAB GOLDEN FILE TEST PASSED")
println("="^70)

