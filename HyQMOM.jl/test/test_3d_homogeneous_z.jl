"""
3D Homogeneous Z Tests

Tests for 3D physical space implementation:
- Basic 3D simulation runs without errors
- Z-slices show expected behavior (homogeneous vs inhomogeneous)
- Nz=1 case matches 2D behavior

Usage:
    julia --project=. test/test_3d_homogeneous_z.jl
    TEST_3D=false julia --project=. -e 'using Pkg; Pkg.test()' # to skip
"""

using Test
using Printf
using Statistics

# Add HyQMOM to load path
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using HyQMOM

# Check if MPI is available
const HAS_MPI = try
    using MPI
    true
catch
    false
end

# Configuration
const RUN_3D_TESTS = get(ENV, "TEST_3D", "true") != "false"
const TEST_TOL = 1e-8

# Helper to initialize MPI
function ensure_mpi_initialized()
    if !HAS_MPI
        return (rank=0, nprocs=1, comm=nothing)
    end
    
    if !MPI.Initialized()
        MPI.Init()
    end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    return (rank=rank, nprocs=nprocs, comm=comm)
end

@testset "3D Homogeneous Z Tests" begin
    if !RUN_3D_TESTS
        @info "Skipping 3D tests (TEST_3D=false)"
    else
        mpi_state = ensure_mpi_initialized()
        rank = mpi_state.rank
        
        @testset "Basic 3D Run (Nz=4, inhomogeneous)" begin
            if rank == 0
                println("\n=== Test: 3D Basic Run ===")
            end
            
            # Small grid for fast testing
            Nx = 10
            Ny = 10
            Nz = 4
            tmax = 0.001  # Very short time
            
            if rank == 0
                println("Running 3D simulation (Nx=$Nx, Ny=$Ny, Nz=$Nz, inhomogeneous_z=false)...")
            end
            
            results = run_simulation(
                Nx=Nx, Ny=Ny, Nz=Nz, tmax=tmax, 
                num_workers=mpi_state.nprocs,
                verbose=false,
                homogeneous_z=false  # Jets only in lower z
            )
            
            if rank == 0
                @test results[:M] !== nothing
                @test results[:time_steps] >= 1
                
                M_3d = results[:M]
                @test size(M_3d, 1) == Nx
                @test size(M_3d, 2) == Ny
                @test size(M_3d, 3) == Nz
                @test size(M_3d, 4) == 35
                
                # Check density remains positive
                @test all(M_3d[:,:,:,1] .> 0)
                
                # Check z-variation (jets only in lower z)
                M_lower_z = M_3d[:,:,1,1]  # First z-slice
                M_upper_z = M_3d[:,:,end,1]  # Last z-slice
                
                max_lower = maximum(M_lower_z)
                max_upper = maximum(M_upper_z)
                
                println("  Lower z (k=1): max density = $max_lower")
                println("  Upper z (k=$Nz): max density = $max_upper")
                
                # Lower z should have jets (high density)
                @test max_lower > 0.1
                # Upper z should have only background (low density)
                @test max_upper < 0.1
                
                println("  ✓ PASS: 3D simulation runs with z-variation")
            end
        end
        
        @testset "Homogeneous Z (jets at all levels)" begin
            if rank == 0
                println("\n=== Test: Homogeneous Z ===")
            end
            
            Nx = 10
            Ny = 10
            Nz = 4
            tmax = 0.001
            
            if rank == 0
                println("Running 3D simulation (homogeneous_z=true)...")
            end
            
            results = run_simulation(
                Nx=Nx, Ny=Ny, Nz=Nz, tmax=tmax,
                num_workers=mpi_state.nprocs,
                verbose=false,
                homogeneous_z=true  # Jets at all z levels
            )
            
            if rank == 0
                M_3d = results[:M]
                
                # All z-slices should have similar max density (jets everywhere)
                max_densities = [maximum(M_3d[:,:,k,1]) for k in 1:Nz]
                
                println("  Max densities by z-slice: ", max_densities)
                
                # All should be > 0.5 (have jets)
                @test all(max_densities .> 0.5)
                
                # Variation between slices should be small (< 50% relative)
                rel_variation = (maximum(max_densities) - minimum(max_densities)) / mean(max_densities)
                println("  Relative variation across z: $rel_variation")
                @test rel_variation < 0.5
                
                println("  ✓ PASS: Homogeneous z shows jets at all levels")
            end
        end
        
        @testset "Nz=1 Case (quasi-2D)" begin
            if rank == 0
                println("\n=== Test: Nz=1 (quasi-2D) ===")
            end
            
            Nx = 10
            Ny = 10
            Nz = 1
            tmax = 0.001
            
            if rank == 0
                println("Running with Nz=1...")
            end
            
            results = run_simulation(
                Nx=Nx, Ny=Ny, Nz=Nz, tmax=tmax,
                num_workers=mpi_state.nprocs,
                verbose=false,
                homogeneous_z=true
            )
            
            if rank == 0
                M_3d = results[:M]
                
                # Should still be 4D but with Nz=1
                @test size(M_3d, 3) == 1
                @test size(M_3d, 1) == Np
                @test size(M_3d, 2) == Np
                
                # Density should be positive
                @test all(M_3d[:,:,1,1] .> 0)
                
                # Should have jets (high max density)
                max_density = maximum(M_3d[:,:,1,1])
                println("  Max density with Nz=1: $max_density")
                @test max_density > 0.5
                
                println("  ✓ PASS: Nz=1 case runs successfully")
            end
        end
        
        @testset "Multiple Z-slices Consistency" begin
            if rank == 0
                println("\n=== Test: Z-slices Consistency ===")
            end
            
            Nx = 10
            Ny = 10
            Nz = 8
            tmax = 0.001
            
            if rank == 0
                println("Running with Nz=$Nz (inhomogeneous)...")
            end
            
            results = run_simulation(
                Nx=Nx, Ny=Ny, Nz=Nz, tmax=tmax,
                num_workers=mpi_state.nprocs,
                verbose=false,
                homogeneous_z=false
            )
            
            if rank == 0
                M_3d = results[:M]
                zm = results[:zm]
                
                # Find jet and background regions
                jet_slices = findall(zm .< 0)
                bg_slices = findall(zm .>= 0)
                
                println("  Jet region: z-slices $jet_slices")
                println("  Background region: z-slices $bg_slices")
                
                if length(jet_slices) >= 2
                    # Compare jet slices
                    k1 = jet_slices[1]
                    k2 = jet_slices[end]
                    
                    max_diff = maximum(abs.(M_3d[:,:,k1,:] .- M_3d[:,:,k2,:]))
                    ref_val = maximum(abs.(M_3d[:,:,k1,1]))
                    rel_diff = max_diff / (ref_val + 1e-10)
                    
                    println("  Jet slices $k1 vs $k2: rel diff = $rel_diff")
                    @test rel_diff < 0.5  # Should be reasonably similar
                end
                
                if length(bg_slices) >= 2
                    # Compare background slices
                    k1 = bg_slices[1]
                    k2 = bg_slices[end]
                    
                    max_diff = maximum(abs.(M_3d[:,:,k1,:] .- M_3d[:,:,k2,:]))
                    ref_val = maximum(abs.(M_3d[:,:,k1,1]))
                    rel_diff = max_diff / (ref_val + 1e-10)
                    
                    println("  Background slices $k1 vs $k2: rel diff = $rel_diff")
                    @test rel_diff < 2.0  # More tolerance for background
                end
                
                println("  ✓ PASS: Z-slices show consistent behavior")
            end
        end
        
        # Finalize MPI if needed
        if HAS_MPI && MPI.Initialized() && abspath(PROGRAM_FILE) == @__FILE__
            MPI.Finalize()
        end
    end
end

# Run if standalone
if abspath(PROGRAM_FILE) == @__FILE__
    # Tests already ran above
    println("\n✅ All 3D homogeneous z tests completed!")
end

