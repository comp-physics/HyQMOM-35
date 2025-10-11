"""
Golden File Validation Tests

These tests compare Julia implementation outputs against MATLAB golden files
to ensure numerical accuracy and correctness of the port.
"""

using Test
using MAT

# Add HyQMOM to load path
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using HyQMOM

const GOLDEN_TOL = 1e-10  # Tolerance for golden file comparisons
const GOLDEN_DIR = joinpath(@__DIR__, "..", "..", "tests", "goldenfiles")

@testset "Golden File Validation" begin
    
    @testset "Moment Conversions Golden" begin
        golden_file = joinpath(GOLDEN_DIR, "test_moment_conversions_golden.mat")
        
        if isfile(golden_file)
            data = matread(golden_file)
            
            # Test 1: Gaussian case
            if haskey(data, "golden_data") && haskey(data["golden_data"], "tests")
                tests = data["golden_data"]["tests"]
                
                if haskey(tests, "gaussian")
                    gaussian = tests["gaussian"]
                    input = gaussian["input"]
                    expected = gaussian["output"]
                    
                    # Run Julia version
                    M = InitializeM4_35(input["rho"], input["u"], input["v"], input["w"],
                                        input["T"], 0.0, 0.0, input["T"], 0.0, input["T"])
                    C4, S4 = M2CS4_35(M)
                    
                    # Compare
                    @test M ≈ vec(expected["M"]) atol=GOLDEN_TOL
                    @test C4 ≈ vec(expected["C4"]) atol=GOLDEN_TOL
                    @test S4 ≈ vec(expected["S4"]) atol=GOLDEN_TOL
                    
                    println("  OK Gaussian case matches golden file")
                end
                
                # Test 2: Correlated case
                if haskey(tests, "correlated")
                    correlated = tests["correlated"]
                    input = correlated["input"]
                    expected = correlated["output"]
                    
                    M = InitializeM4_35(input["rho"], input["u"], input["v"], input["w"],
                                        input["C200"], input["C110"], input["C101"],
                                        input["C020"], input["C011"], input["C002"])
                    C4, S4 = M2CS4_35(M)
                    
                    @test M ≈ vec(expected["M"]) atol=GOLDEN_TOL
                    @test C4 ≈ vec(expected["C4"]) atol=GOLDEN_TOL
                    @test S4 ≈ vec(expected["S4"]) atol=GOLDEN_TOL
                    
                    println("  OK Correlated case matches golden file")
                end
            end
        else
            @warn "Golden file not found: $golden_file"
        end
    end
    
    @testset "Initialization Golden" begin
        golden_file = joinpath(GOLDEN_DIR, "test_initialization_golden.mat")
        
        if isfile(golden_file)
            data = matread(golden_file)
            
            if haskey(data, "golden_data") && haskey(data["golden_data"], "tests")
                tests = data["golden_data"]["tests"]
                
                # Test all cases
                for (case_name, case_data) in tests
                    if haskey(case_data, "input") && haskey(case_data, "output")
                        input = case_data["input"]
                        expected = vec(case_data["output"])
                        
                        # Run Julia version
                        if all(input["corr"] .== 0)
                            M = InitializeM4_35(input["rho"], input["u"], input["v"], input["w"],
                                                input["T"], 0.0, 0.0, input["T"], 0.0, input["T"])
                        else
                            T = input["T"]
                            C200, C020, C002 = T, T, T
                            C110 = input["corr"][1] * sqrt(C200 * C020)
                            C101 = input["corr"][2] * sqrt(C200 * C002)
                            C011 = input["corr"][3] * sqrt(C020 * C002)
                            M = InitializeM4_35(input["rho"], input["u"], input["v"], input["w"],
                                                C200, C110, C101, C020, C011, C002)
                        end
                        
                        @test M ≈ expected atol=GOLDEN_TOL
                        println("  OK $case_name matches golden file")
                    end
                end
            end
        else
            @warn "Golden file not found: $golden_file"
        end
    end
    
    @testset "Realizability Golden" begin
        golden_file = joinpath(GOLDEN_DIR, "test_realizability_golden.mat")
        
        if isfile(golden_file)
            data = matread(golden_file)
            
            if haskey(data, "golden_data") && haskey(data["golden_data"], "tests")
                tests = data["golden_data"]["tests"]
                
                # Test S2 realizability cases
                for (case_name, case_data) in tests
                    if startswith(case_name, "S2_case") && haskey(case_data, "input")
                        input = case_data["input"]
                        expected = case_data["output"]
                        
                        # Run Julia version
                        S110r, S101r, S011r, S2r = realizability(:S2, input["S110"], 
                                                                   input["S101"], input["S011"])
                        
                        @test S110r ≈ expected["S110r"] atol=GOLDEN_TOL
                        @test S101r ≈ expected["S101r"] atol=GOLDEN_TOL
                        @test S011r ≈ expected["S011r"] atol=GOLDEN_TOL
                        @test S2r ≈ expected["S2r"] atol=GOLDEN_TOL
                        
                        println("  OK $case_name matches golden file")
                    end
                end
            end
        else
            @warn "Golden file not found: $golden_file"
        end
    end
    
    @testset "Closures Golden" begin
        golden_file = joinpath(GOLDEN_DIR, "test_closures_golden.mat")
        
        if isfile(golden_file)
            data = matread(golden_file)
            
            if haskey(data, "golden_data") && haskey(data["golden_data"], "tests")
                tests = data["golden_data"]["tests"]
                
                # Test Gaussian case
                if haskey(tests, "gaussian")
                    gaussian = tests["gaussian"]
                    input = vec(gaussian["input"])
                    expected = gaussian["output"]
                    
                    # Run Julia version
                    M5, C5, S5 = Moments5_3D(input)
                    
                    @test M5 ≈ vec(expected["M5"]) atol=GOLDEN_TOL
                    @test C5 ≈ vec(expected["C5"]) atol=GOLDEN_TOL
                    @test S5 ≈ vec(expected["S5"]) atol=GOLDEN_TOL
                    
                    println("  OK Gaussian Moments5_3D matches golden file")
                end
                
                # Test 1D closure
                if haskey(tests, "closure_1d")
                    closure_1d = tests["closure_1d"]
                    input = vec(closure_1d["input"])
                    expected = closure_1d["output"]
                    
                    # Run Julia version
                    Mp, vpmin, vpmax = closure_and_eigenvalues(input)
                    
                    # Mp is a scalar, not a vector
                    expected_Mp = expected["Mp"]
                    if expected_Mp isa Array
                        expected_Mp = Float64(expected_Mp[1])
                    else
                        expected_Mp = Float64(expected_Mp)
                    end
                    
                    @test Mp ≈ expected_Mp atol=GOLDEN_TOL
                    @test vpmin ≈ Float64(expected["vpmin"]) atol=GOLDEN_TOL
                    @test vpmax ≈ Float64(expected["vpmax"]) atol=GOLDEN_TOL
                    
                    println("  OK 1D closure matches golden file")
                end
            end
        else
            @warn "Golden file not found: $golden_file"
        end
    end
    
    @testset "Flux and Eigenvalues Golden" begin
        golden_file = joinpath(GOLDEN_DIR, "test_flux_eigenvalues_golden.mat")
        
        if isfile(golden_file)
            data = matread(golden_file)
            
            if haskey(data, "golden_data") && haskey(data["golden_data"], "tests")
                tests = data["golden_data"]["tests"]
                
                # Test flux computation
                if haskey(tests, "flux_basic")
                    flux_basic = tests["flux_basic"]
                    input = vec(flux_basic["input"])
                    expected = flux_basic["output"]
                    
                    # Run Julia version
                    Fx, Fy, Fz, M_real = Flux_closure35_and_realizable_3D(input, 1, 0.5)
                    
                    @test Fx ≈ vec(expected["Fx"]) atol=GOLDEN_TOL
                    @test Fy ≈ vec(expected["Fy"]) atol=GOLDEN_TOL
                    @test Fz ≈ vec(expected["Fz"]) atol=GOLDEN_TOL
                    @test M_real ≈ vec(expected["M_real"]) atol=GOLDEN_TOL
                    
                    println("  OK Flux computation matches golden file")
                end
                
                # Test eigenvalues
                if haskey(tests, "eigenvalues")
                    eigenvalues_test = tests["eigenvalues"]
                    input = vec(eigenvalues_test["input"])
                    expected = eigenvalues_test["output"]
                    
                    # Run Julia version
                    v6min_x, v6max_x = eigenvalues6_hyperbolic_3D(input, 1, 1, 0.5)
                    v6min_y, v6max_y = eigenvalues6_hyperbolic_3D(input, 2, 1, 0.5)
                    
                    @test v6min_x ≈ Float64(expected["v6min_x"]) atol=GOLDEN_TOL
                    @test v6max_x ≈ Float64(expected["v6max_x"]) atol=GOLDEN_TOL
                    @test v6min_y ≈ Float64(expected["v6min_y"]) atol=GOLDEN_TOL
                    @test v6max_y ≈ Float64(expected["v6max_y"]) atol=GOLDEN_TOL
                    
                    println("  OK Eigenvalues match golden file")
                end
            end
        else
            @warn "Golden file not found: $golden_file"
        end
    end
    
    @testset "Numerical Schemes Golden" begin
        golden_file = joinpath(GOLDEN_DIR, "test_numerical_schemes_golden.mat")
        
        if isfile(golden_file)
            data = matread(golden_file)
            
            if haskey(data, "golden_data") && haskey(data["golden_data"], "tests")
                tests = data["golden_data"]["tests"]
                
                # Test pas_HLL
                if haskey(tests, "pas_hll")
                    pas_hll_test = tests["pas_hll"]
                    input = pas_hll_test["input"]
                    expected = pas_hll_test["output"]
                    
                    # Run Julia version
                    Mp = pas_HLL(input["M"], input["F"], input["dt"], input["dx"],
                                 vec(input["vpmin"]), vec(input["vpmax"]);
                                 apply_bc_left=true, apply_bc_right=true)
                    
                    @test Mp ≈ expected atol=GOLDEN_TOL
                    
                    println("  OK pas_HLL matches golden file")
                end
            end
        else
            @warn "Golden file not found: $golden_file"
        end
    end
    
    @testset "Collision Golden" begin
        golden_file = joinpath(GOLDEN_DIR, "test_collision_golden.mat")
        
        if isfile(golden_file)
            data = matread(golden_file)
            
            if haskey(data, "golden_data") && haskey(data["golden_data"], "tests")
                tests = data["golden_data"]["tests"]
                
                # Test collision35
                if haskey(tests, "basic")
                    basic = tests["basic"]
                    input = basic["input"]
                    expected = vec(basic["output"])
                    
                    # Run Julia version
                    M_out = collision35(vec(input["M"]), input["dt"], input["Kn"])
                    
                    @test M_out ≈ expected atol=GOLDEN_TOL
                    
                    println("  OK collision35 matches golden file")
                end
            end
        else
            @warn "Golden file not found: $golden_file"
        end
    end
    
end
