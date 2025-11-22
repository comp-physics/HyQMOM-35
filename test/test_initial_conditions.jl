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
        # Create simple grid parameters
        Nx, Ny, Nz = 10, 10, 1
        xmin, xmax = -0.5, 0.5
        ymin, ymax = -0.5, 0.5
        zmin, zmax = -0.5, 0.5
        
        xm = collect(range(xmin, xmax, length=Nx))
        ym = collect(range(ymin, ymax, length=Ny))
        zm = collect(range(zmin, zmax, length=Nz))
        
        grid_params = (Nx=Nx, Ny=Ny, Nz=Nz, 
                       xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                       xm=xm, ym=ym, zm=zm)
        
        # Background
        background = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (Inf, Inf, Inf),
            density = 0.1,
            velocity = (0.0, 0.0, 0.0),
            temperature = 1.0
        )
        
        # Single jet
        jet = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (0.2, 0.2, 0.2),
            density = 1.0,
            velocity = (0.5, 0.0, 0.0),
            temperature = 1.0
        )
        
        M = HyQMOM.initialize_moment_field(grid_params, background, [jet])
        
        # Check dimensions
        @test size(M) == (Nx, Ny, Nz, 35)
        
        # Check all values are finite
        @test all(isfinite, M)
        
        # Check that density is positive everywhere
        @test all(M[:, :, :, 1] .> 0)
        
        # Check that jet region has higher density than background
        center_idx = div(Nx, 2) + 1
        @test M[center_idx, center_idx, 1, 1] > 0.5  # Should be closer to jet density
        
        # Test with no regions (just background)
        M_bg = HyQMOM.initialize_moment_field(grid_params, background, HyQMOM.CubicRegion[])
        @test size(M_bg) == (Nx, Ny, Nz, 35)
        @test all(isfinite, M_bg)
        @test all(M_bg[:, :, :, 1] .≈ 0.1)  # All should be background density
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
    
    @testset "place_cubic_region!" begin
        Nx, Ny, Nz = 10, 10, 1
        xmin, xmax = -0.5, 0.5
        ymin, ymax = -0.5, 0.5
        zmin, zmax = -0.5, 0.5
        
        xm = collect(range(xmin, xmax, length=Nx))
        ym = collect(range(ymin, ymax, length=Ny))
        zm = collect(range(zmin, zmax, length=Nz))
        
        M = zeros(Float64, Nx, Ny, Nz, 35)
        
        # Place a region
        region = HyQMOM.CubicRegion(
            center = (0.0, 0.0, 0.0),
            width = (0.2, 0.2, 0.2),
            density = 1.0,
            velocity = (0.5, 0.0, 0.0),
            temperature = 1.0
        )
        
        HyQMOM.place_cubic_region!(M, region, xm, ym, zm, 0.0, 0.0, 0.0)
        
        # Check that some cells were modified
        @test any(M[:, :, :, 1] .> 0)
        
        # Check center cell should have jet density
        center_idx = div(Nx, 2) + 1
        @test M[center_idx, center_idx, 1, 1] ≈ 1.0
        
        # Test with use_overlap=false
        M2 = zeros(Float64, Nx, Ny, Nz, 35)
        HyQMOM.place_cubic_region!(M2, region, xm, ym, zm, 0.0, 0.0, 0.0, use_overlap=false)
        @test any(M2[:, :, :, 1] .> 0)
    end
    
    @testset "Full initialization workflow" begin
        # Test complete workflow with crossing jets
        Nx, Ny, Nz = 20, 20, 1
        xmin, xmax = -0.5, 0.5
        ymin, ymax = -0.5, 0.5
        zmin, zmax = -0.5, 0.5
        
        xm = collect(range(xmin, xmax, length=Nx))
        ym = collect(range(ymin, ymax, length=Ny))
        zm = collect(range(zmin, zmax, length=Nz))
        
        grid_params = (Nx=Nx, Ny=Ny, Nz=Nz, 
                       xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                       xm=xm, ym=ym, zm=zm)
        
        for Ma in [0.0, 1.0, 5.0, 10.0]
            for Kn in [0.1, 1.0]
                background, jets = HyQMOM.crossing_jets_ic(
                    Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax,
                    Ma=Ma
                )
                
                M = HyQMOM.initialize_moment_field(grid_params, background, jets)
                
                # Basic sanity checks
                @test size(M) == (Nx, Ny, Nz, 35)
                @test all(isfinite, M)
                @test all(M[:, :, :, 1] .> 0)  # Positive density
                
                # Check that we have variation (not all same)
                @test maximum(M[:, :, :, 1]) > minimum(M[:, :, :, 1])
            end
        end
    end
    
end

