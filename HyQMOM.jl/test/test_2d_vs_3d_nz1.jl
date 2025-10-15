"""
Test: Validate 3D code with Nz=1 matches archived 2D code

This test runs both the archived 2D Julia code and the new 3D Julia code
with Nz=1, and compares the results to ensure backward compatibility.
"""

using Test
using Printf

# Test configuration
const Nx = 20
const Ny = 20
const tmax = 0.01  # Short simulation for quick comparison
const TEST_TOL_ABS = 1e-10
const TEST_TOL_REL = 1e-8

println("="^70)
println("Testing 3D (Nz=1) vs Archived 2D Code")
println("="^70)

# Add paths
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "HyQMOM_2D_archive"))

# Load 3D version (current)
println("\nLoading 3D HyQMOM module...")
import HyQMOM
const HyQMOM_3D = HyQMOM

# Load 2D version (archived)
println("Loading 2D HyQMOM module (archived)...")
module HyQMOM_2D_Wrapper
    include("../../HyQMOM_2D_archive/src/HyQMOM.jl")
end
const HyQMOM_2D = HyQMOM_2D_Wrapper.HyQMOM

println("✓ Both modules loaded successfully")

@testset "2D vs 3D (Nz=1) Comparison" begin
    
    @testset "Setup and Run Simulations" begin
        println("\n" * "="^70)
        println("Running 2D simulation (archived code)...")
        println("="^70)
        
        # Run 2D simulation
        params_2d = (
            Nx=Nx, Ny=Ny, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
            dx=1.0/Np, dy=1.0/Np, Nmom=35, nnmax=20000000, dtmax=1.0,
            rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
            symmetry_check_interval=10,
            enable_memory_tracking=false,
            debug_output=false
        )
        
        M_2d, final_time_2d, time_steps_2d, grid_2d = HyQMOM_2D.simulation_runner(params_2d)
        
        @test M_2d !== nothing
        @test size(M_2d, 1) == Np
        @test size(M_2d, 2) == Np
        @test size(M_2d, 3) == 35
        @test ndims(M_2d) == 3
        
        println("\n2D Results:")
        println("  Final time: $final_time_2d")
        println("  Time steps: $time_steps_2d")
        println("  M shape: $(size(M_2d))")
        
        println("\n" * "="^70)
        println("Running 3D simulation with Nz=1 (current code)...")
        println("="^70)
        
        # Run 3D simulation with Nz=1
        params_3d = (
            Nx=Nx, Ny=Ny, Nz=1, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
            dx=1.0/Np, dy=1.0/Np, dz=1.0, Nmom=35, nnmax=20000000, dtmax=1.0,
            rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
            symmetry_check_interval=10,
            homogeneous_z=true,  # Irrelevant for Nz=1
            enable_memory_tracking=false,
            debug_output=false
        )
        
        M_3d, final_time_3d, time_steps_3d, grid_3d = HyQMOM_3D.simulation_runner(params_3d)
        
        @test M_3d !== nothing
        @test size(M_3d, 1) == Np
        @test size(M_3d, 2) == Np
        @test size(M_3d, 3) == 1  # Nz=1
        @test size(M_3d, 4) == 35
        @test ndims(M_3d) == 4
        
        println("\n3D Results:")
        println("  Final time: $final_time_3d")
        println("  Time steps: $time_steps_3d")
        println("  M shape: $(size(M_3d))")
        
        # Extract the single z-slice from 3D
        M_3d_slice = M_3d[:, :, 1, :]
        
        @test size(M_3d_slice) == size(M_2d)
        println("\n✓ Array dimensions match after extracting z-slice")
    end
    
    @testset "Compare Final Times" begin
        println("\n" * "="^70)
        println("Comparing final times...")
        println("="^70)
        
        # Re-run to get fresh results (variables from nested testset)
        params_2d = (Np=Np, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
                    dx=1.0/Np, dy=1.0/Np, Nmom=35, nnmax=20000000, dtmax=1.0,
                    rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
                    symmetry_check_interval=10, enable_memory_tracking=false, debug_output=false)
        M_2d, final_time_2d, time_steps_2d, grid_2d = HyQMOM_2D.simulation_runner(params_2d)
        
        params_3d = (Np=Np, Nz=1, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
                    dx=1.0/Np, dy=1.0/Np, dz=1.0, Nmom=35, nnmax=20000000, dtmax=1.0,
                    rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
                    symmetry_check_interval=10, homogeneous_z=true,
                    enable_memory_tracking=false, debug_output=false)
        M_3d, final_time_3d, time_steps_3d, grid_3d = HyQMOM_3D.simulation_runner(params_3d)
        
        M_3d_slice = M_3d[:, :, 1, :]
        
        time_diff = abs(final_time_3d - final_time_2d)
        println("  2D final time: $final_time_2d")
        println("  3D final time: $final_time_3d")
        println("  Difference: $time_diff")
        
        @test abs(final_time_3d - final_time_2d) < 1e-10
        
        if time_steps_2d == time_steps_3d
            println("  ✓ Time step counts match: $time_steps_2d")
        else
            @warn "Time step counts differ: 2D=$time_steps_2d, 3D=$time_steps_3d"
            # This is OK if they both reach tmax - may use slightly different dt
        end
    end
    
    @testset "Compare Moment Arrays" begin
        println("\n" * "="^70)
        println("Comparing moment arrays...")
        println("="^70)
        
        # Re-run to get results
        params_2d = (Np=Np, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
                    dx=1.0/Np, dy=1.0/Np, Nmom=35, nnmax=20000000, dtmax=1.0,
                    rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
                    symmetry_check_interval=10, enable_memory_tracking=false, debug_output=false)
        M_2d, final_time_2d, time_steps_2d, grid_2d = HyQMOM_2D.simulation_runner(params_2d)
        
        params_3d = (Np=Np, Nz=1, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
                    dx=1.0/Np, dy=1.0/Np, dz=1.0, Nmom=35, nnmax=20000000, dtmax=1.0,
                    rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
                    symmetry_check_interval=10, homogeneous_z=true,
                    enable_memory_tracking=false, debug_output=false)
        M_3d, final_time_3d, time_steps_3d, grid_3d = HyQMOM_3D.simulation_runner(params_3d)
        
        M_3d_slice = M_3d[:, :, 1, :]
        
        # Compute differences
        abs_diff = abs.(M_3d_slice .- M_2d)
        rel_diff = abs_diff ./ (abs.(M_2d) .+ 1e-30)
        
        max_abs_diff = maximum(abs_diff)
        max_rel_diff = maximum(rel_diff)
        
        # Find location of max difference
        idx_max = argmax(abs_diff)
        i_max, j_max, m_max = Tuple(idx_max)
        
        println("\nDifference Statistics:")
        println("  Max absolute difference: $max_abs_diff")
        println("  Max relative difference: $max_rel_diff")
        println("  Location: (i=$i_max, j=$j_max, moment=$m_max)")
        println("  2D value: $(M_2d[i_max, j_max, m_max])")
        println("  3D value: $(M_3d_slice[i_max, j_max, m_max])")
        
        # Check key moments
        println("\nChecking key moments:")
        for (m, name) in [(1, "M000"), (2, "M100"), (6, "M010"), (3, "M200"), (10, "M020")]
            diff_m = maximum(abs.(M_3d_slice[:, :, m] .- M_2d[:, :, m]))
            println("  $name: max diff = $diff_m")
        end
        
        # Test tolerances
        @test max_abs_diff < TEST_TOL_ABS || max_rel_diff < TEST_TOL_REL
        
        if max_abs_diff < TEST_TOL_ABS && max_rel_diff < TEST_TOL_REL
            println("\n✓ 3D (Nz=1) matches 2D within tolerance!")
            println("  Absolute: $max_abs_diff < $TEST_TOL_ABS")
            println("  Relative: $max_rel_diff < $TEST_TOL_REL")
        elseif max_abs_diff < TEST_TOL_ABS
            println("\n✓ 3D (Nz=1) matches 2D within absolute tolerance")
            println("  Absolute: $max_abs_diff < $TEST_TOL_ABS")
            @warn "Relative difference large but absolute is OK: $max_rel_diff"
        else
            println("\n✗ Differences exceed tolerance!")
            @warn "This may indicate a problem with the 3D implementation"
        end
    end
    
    @testset "Compare Grid Structures" begin
        println("\n" * "="^70)
        println("Comparing grid structures...")
        println("="^70)
        
        # Re-run simulations
        params_2d = (Np=Np, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
                    dx=1.0/Np, dy=1.0/Np, Nmom=35, nnmax=20000000, dtmax=1.0,
                    rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
                    symmetry_check_interval=10, enable_memory_tracking=false, debug_output=false)
        M_2d, final_time_2d, time_steps_2d, grid_2d = HyQMOM_2D.simulation_runner(params_2d)
        
        params_3d = (Np=Np, Nz=1, tmax=tmax, Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5,
                    dx=1.0/Np, dy=1.0/Np, dz=1.0, Nmom=35, nnmax=20000000, dtmax=1.0,
                    rhol=1.0, rhor=0.01, T=1.0, r110=0.0, r101=0.0, r011=0.0,
                    symmetry_check_interval=10, homogeneous_z=true,
                    enable_memory_tracking=false, debug_output=false)
        M_3d, final_time_3d, time_steps_3d, grid_3d = HyQMOM_3D.simulation_runner(params_3d)
        
        # Compare x and y grids
        @test length(grid_2d.xm) == length(grid_3d.xm)
        @test length(grid_2d.ym) == length(grid_3d.ym)
        @test maximum(abs.(grid_2d.xm .- grid_3d.xm)) < 1e-15
        @test maximum(abs.(grid_2d.ym .- grid_3d.ym)) < 1e-15
        
        # Check z grid
        @test haskey(grid_3d, :z)
        @test haskey(grid_3d, :zm)
        @test length(grid_3d.z) == 2  # Cell edges: 2 points for 1 cell
        @test length(grid_3d.zm) == 1  # Cell centers: 1 point for 1 cell
        
        println("  ✓ X and Y grids match")
        println("  ✓ Z grid properly defined for Nz=1")
    end
end

println("\n" * "="^70)
println("Test Complete!")
println("="^70)

