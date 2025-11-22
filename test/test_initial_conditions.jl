"""
Tests for initial_conditions.jl

Tests the flexible initial condition system including:
- CubicRegion construction
- Moment field initialization
- Region placement and overlap detection
- Crossing jets convenience function
"""

using Test
using HyQMOM

@testset "Initial Conditions" begin
    
    @testset "CubicRegion construction" begin
        # Test basic construction
        region = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (0.2, 0.2, 0.2),
            density = 1.0,
            velocity = (0.5, 0.0, 0.0),
            temperature = 1.0
        )
        
        @test region.center == (0.0, 0.0, 0.0)
        @test region.width == (0.2, 0.2, 0.2)
        @test region.density == 1.0
        @test region.velocity == (0.5, 0.0, 0.0)
        @test region.temperature == 1.0
        
        # Test with different inputs (should convert to Float64)
        region2 = HyQMOM.CubicRegion(
            center = [1, 2, 3],
            width = [0.1, 0.1, 0.1],
            density = 2,
            velocity = [0, 0, 1],
            temperature = 2
        )
        
        @test region2.center isa NTuple{3, Float64}
        @test region2.width isa NTuple{3, Float64}
        @test region2.density isa Float64
        @test region2.velocity isa NTuple{3, Float64}
        @test region2.temperature isa Float64
        
        # Test default temperature
        region3 = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (0.1, 0.1, 0.1),
            density = 1.0,
            velocity = (0.0, 0.0, 0.0)
        )
        @test region3.temperature == 1.0
    end
    
    @testset "region_to_moments" begin
        # Test stationary region with unit temperature
        region = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (1.0, 1.0, 1.0),
            density = 1.0,
            velocity = (0.0, 0.0, 0.0),
            temperature = 1.0
        )
        
        Mr = HyQMOM.region_to_moments(region, 0.0, 0.0, 0.0)
        
        @test length(Mr) == 35
        @test all(isfinite, Mr)
        @test Mr[1] ≈ 1.0  # Density should be 1.0
        
        # Test moving region
        region_moving = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (1.0, 1.0, 1.0),
            density = 2.0,
            velocity = (1.0, 0.5, 0.0),
            temperature = 0.5
        )
        
        Mr_moving = HyQMOM.region_to_moments(region_moving, 0.0, 0.0, 0.0)
        
        @test length(Mr_moving) == 35
        @test all(isfinite, Mr_moving)
        @test Mr_moving[1] ≈ 2.0  # Density
        
        # Test with correlations
        Mr_corr = HyQMOM.region_to_moments(region, 0.1, 0.2, 0.3)
        @test length(Mr_corr) == 35
        @test all(isfinite, Mr_corr)
    end
    
    @testset "point_in_cube" begin
        cube = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (1.0, 1.0, 1.0),
            density = 1.0,
            velocity = (0.0, 0.0, 0.0),
            temperature = 1.0
        )
        
        # Points inside
        @test HyQMOM.point_in_cube((0.0, 0.0, 0.0), cube) == true
        @test HyQMOM.point_in_cube((0.25, 0.25, 0.25), cube) == true
        @test HyQMOM.point_in_cube((-0.49, 0.0, 0.0), cube) == true
        @test HyQMOM.point_in_cube((0.0, 0.49, 0.0), cube) == true
        @test HyQMOM.point_in_cube((0.0, 0.0, 0.49), cube) == true
        
        # Points outside
        @test HyQMOM.point_in_cube((0.51, 0.0, 0.0), cube) == false
        @test HyQMOM.point_in_cube((0.0, 0.51, 0.0), cube) == false
        @test HyQMOM.point_in_cube((0.0, 0.0, 0.51), cube) == false
        @test HyQMOM.point_in_cube((1.0, 1.0, 1.0), cube) == false
        
        # Test infinite width (background)
        background = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (Inf, Inf, Inf),
            density = 0.01,
            velocity = (0.0, 0.0, 0.0),
            temperature = 1.0
        )
        
        @test HyQMOM.point_in_cube((100.0, 100.0, 100.0), background) == true
        @test HyQMOM.point_in_cube((-100.0, -100.0, -100.0), background) == true
    end
    
    @testset "cell_overlaps_cube" begin
        cube = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (0.2, 0.2, 0.2),
            density = 1.0,
            velocity = (0.0, 0.0, 0.0),
            temperature = 1.0
        )
        
        cell_size = (0.1, 0.1, 0.1)
        
        # Cell overlaps cube
        @test HyQMOM.cell_overlaps_cube((0.0, 0.0, 0.0), cell_size, cube) == true
        @test HyQMOM.cell_overlaps_cube((0.05, 0.05, 0.05), cell_size, cube) == true
        
        # Cell just touching cube (should overlap)
        @test HyQMOM.cell_overlaps_cube((0.15, 0.0, 0.0), cell_size, cube) == true
        
        # Cell completely outside
        @test HyQMOM.cell_overlaps_cube((0.5, 0.0, 0.0), cell_size, cube) == false
        @test HyQMOM.cell_overlaps_cube((0.0, 0.5, 0.0), cell_size, cube) == false
        
        # Test with infinite width
        background = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (Inf, Inf, Inf),
            density = 0.01,
            velocity = (0.0, 0.0, 0.0),
            temperature = 1.0
        )
        
        @test HyQMOM.cell_overlaps_cube((100.0, 100.0, 100.0), cell_size, background) == true
    end
    
    @testset "initialize_moment_field" begin
        # Test basic initialization
        Nx, Ny, Nz = 8, 8, 4
        xmin, xmax = -1.0, 1.0
        ymin, ymax = -1.0, 1.0
        zmin, zmax = -0.5, 0.5
        
        xm = range(xmin, xmax, length=Nx)
        ym = range(ymin, ymax, length=Ny)
        zm = range(zmin, zmax, length=Nz)
        
        grid_params = (Nx=Nx, Ny=Ny, Nz=Nz, xm=xm, ym=ym, zm=zm,
                      xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
        
        background = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (Inf, Inf, Inf),
            density = 0.1,
            velocity = (0.0, 0.0, 0.0),
            temperature = 1.0
        )
        
        region1 = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (0.5, 0.5, 0.5),
            density = 1.0,
            velocity = (0.5, 0.0, 0.0),
            temperature = 1.5
        )
        
        M = HyQMOM.initialize_moment_field(grid_params, background, [region1])
        
        @test size(M) == (Nx, Ny, Nz, 35)
        @test all(isfinite, M)
        
        # Check that background density appears in corners
        @test M[1,1,1,1] ≈ 0.1 || M[1,1,1,1] ≈ 1.0  # Could be either background or region
        
        # Test with multiple regions
        region2 = HyQMOM.CubicRegion(
            center = (0.5, 0.5, 0.0),
            width = (0.3, 0.3, 0.3),
            density = 2.0,
            velocity = (-0.5, -0.5, 0.0),
            temperature = 0.8
        )
        
        M2 = HyQMOM.initialize_moment_field(grid_params, background, [region1, region2])
        @test size(M2) == (Nx, Ny, Nz, 35)
        @test all(isfinite, M2)
        
        # Test with correlations
        M3 = HyQMOM.initialize_moment_field(grid_params, background, [region1], 
                                           r110=0.1, r101=0.2, r011=0.3)
        @test size(M3) == (Nx, Ny, Nz, 35)
        @test all(isfinite, M3)
    end
    
    @testset "place_cubic_region!" begin
        Nx, Ny, Nz = 10, 10, 5
        xm = range(-1.0, 1.0, length=Nx)
        ym = range(-1.0, 1.0, length=Ny)
        zm = range(-0.5, 0.5, length=Nz)
        
        M = zeros(Nx, Ny, Nz, 35)
        
        # Fill with background
        bg_moments = HyQMOM.region_to_moments(
            HyQMOM.CubicRegion(center=(0.0,0.0,0.0), width=(Inf,Inf,Inf), 
                             density=0.1, velocity=(0.0,0.0,0.0), temperature=1.0),
            0.0, 0.0, 0.0
        )
        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            M[i,j,k,:] = bg_moments
        end
        
        # Place a region with overlap detection (default)
        region = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (0.4, 0.4, 0.4),
            density = 1.0,
            velocity = (0.5, 0.0, 0.0),
            temperature = 1.0
        )
        
        HyQMOM.place_cubic_region!(M, region, xm, ym, zm, 0.0, 0.0, 0.0, use_overlap=true)
        
        @test any(M[:,:,:,1] .≈ 1.0)  # Should have placed some cells with density 1.0
        @test any(M[:,:,:,1] .≈ 0.1)  # Should still have background in some cells
        
        # Test legacy point-in-cube mode
        M_legacy = zeros(Nx, Ny, Nz, 35)
        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            M_legacy[i,j,k,:] = bg_moments
        end
        
        HyQMOM.place_cubic_region!(M_legacy, region, xm, ym, zm, 0.0, 0.0, 0.0, use_overlap=false)
        
        @test any(M_legacy[:,:,:,1] .≈ 1.0)  # Should have placed some cells
        @test all(isfinite, M_legacy)
    end
    
    @testset "crossing_jets_ic" begin
        Nx, Ny, Nz = 20, 20, 1
        xmin, xmax = -0.5, 0.5
        ymin, ymax = -0.5, 0.5
        zmin, zmax = -0.5, 0.5
        
        # Test Ma=0 case
        background, jets = HyQMOM.crossing_jets_ic(
            Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax,
            Ma=0.0, rhol=1.0, rhor=0.01, T=1.0
        )
        
        @test background.density ≈ 0.01
        @test background.velocity == (0.0, 0.0, 0.0)
        @test length(jets) == 2
        @test jets[1].density ≈ 1.0
        @test jets[2].density ≈ 1.0
        @test jets[1].velocity == (0.0, 0.0, 0.0)  # Ma=0 means stationary
        @test jets[2].velocity == (0.0, 0.0, 0.0)
        
        # Test Ma=5 case
        background2, jets2 = HyQMOM.crossing_jets_ic(
            Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax,
            Ma=5.0
        )
        
        Uc_expected = 5.0 / sqrt(2.0)
        @test jets2[1].velocity[1] ≈ Uc_expected
        @test jets2[1].velocity[2] ≈ Uc_expected
        @test jets2[2].velocity[1] ≈ -Uc_expected
        @test jets2[2].velocity[2] ≈ -Uc_expected
        
        # Test custom jet size
        background3, jets3 = HyQMOM.crossing_jets_ic(
            Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax,
            jet_size=0.2
        )
        
        jet_width_expected = 0.2 * (xmax - xmin)
        @test jets3[1].width[1] ≈ jet_width_expected
        @test jets3[1].width[2] ≈ jet_width_expected
        
        # Test jet offset
        background4, jets4 = HyQMOM.crossing_jets_ic(
            Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax,
            jet_size=0.1, jet_offset=1.0
        )
        
        # Jets should be offset from center
        @test jets4[1].center[1] != 0.0
        @test jets4[1].center[2] != 0.0
    end
end

