"""
Unit tests for initial condition properties.

Tests the geometric accuracy and physical properties of initial conditions.
"""

using Test
using HyQMOM

@testset "Initial Condition Properties" begin
    
    @testset "CubicRegion basic properties" begin
        # Test CubicRegion constructor
        center = (0.5, 0.5, 0.5)
        width = (0.2, 0.2, 0.2)
        density = 2.0
        velocity = (10.0, 5.0, -3.0)
        temperature = 1.0
        
        jet = HyQMOM.CubicRegion(
            center=center,
            width=width,
            density=density,
            velocity=velocity,
            temperature=temperature
        )
        
        @test jet.center == center
        @test jet.width == width
        @test jet.density == density
        @test jet.velocity == velocity
        @test jet.temperature == temperature
    end
    
    @testset "point_in_cube basic test" begin
        # Test point_in_cube function with simple cases
        center = (0.5, 0.5, 0.5)
        width = (0.2, 0.2, 0.2)  # half_width = 0.1 each direction
        
        jet = HyQMOM.CubicRegion(
            center=center,
            width=width,
            density=2.0,
            velocity=(0.0, 0.0, 0.0),
            temperature=1.0
        )
        
        # Point at center should be inside
        @test HyQMOM.point_in_cube((0.5, 0.5, 0.5), jet)
        
        # Point at corner should be inside
        @test HyQMOM.point_in_cube((0.6, 0.6, 0.6), jet)
        
        # Point far outside should not be inside
        @test !HyQMOM.point_in_cube((0.8, 0.8, 0.8), jet)
    end
    
    @testset "point_in_cube boundary test" begin
        # Test exact boundaries
        center = (0.5, 0.5, 0.5)
        width = (0.2, 0.2, 0.2)
        
        jet = HyQMOM.CubicRegion(
            center=center,
            width=width,
            density=2.0,
            velocity=(0.0, 0.0, 0.0),
            temperature=1.0
        )
        
        # Point exactly at edge (x_max) - width/2 = 0.1 from center
        edge_point = (0.6, 0.5, 0.5)
        @test HyQMOM.point_in_cube(edge_point, jet)
        
        # Point just outside edge
        outside_point = (0.61, 0.5, 0.5)
        @test !HyQMOM.point_in_cube(outside_point, jet)
    end
    
    @testset "Sharp corners (no rounding)" begin
        # Verify corners are sharp, not rounded
        center = (0.5, 0.5, 0.5)
        width = (0.2, 0.2, 0.2)
        
        jet = HyQMOM.CubicRegion(
            center=center,
            width=width,
            density=2.0,
            velocity=(0.0, 0.0, 0.0),
            temperature=1.0
        )
        
        # Corner point should be inside (sharp corner)
        corner = (0.6, 0.6, 0.6)
        @test HyQMOM.point_in_cube(corner, jet)
        
        # Just outside corner should be outside
        outside_corner = (0.61, 0.61, 0.61)
        @test !HyQMOM.point_in_cube(outside_corner, jet)
    end
    
    @testset "Jet placement accuracy" begin
        # Test that jets are placed at correct location in field
        # This is a conceptual test for the initialize_moment_field logic
        
        Lx, Ly, Lz = 1.0, 1.0, 1.0
        Nx, Ny, Nz = 20, 20, 20
        dx = Lx / Nx
        dy = Ly / Ny
        dz = Lz / Nz
        
        # Jet at center
        center = (0.5, 0.5, 0.5)
        half_width = 0.1
        
        # Cell index at center
        i_center = round(Int, center[1] / dx)
        j_center = round(Int, center[2] / dy)
        k_center = round(Int, center[3] / dz)
        
        # Should be roughly in middle of domain
        @test 5 < i_center < 15
        @test 5 < j_center < 15
        @test 5 < k_center < 15
    end
    
    @testset "Background fills correctly" begin
        # Test that background state is set up correctly
        rho_bg = 0.5
        u_bg, v_bg, w_bg = 0.0, 0.0, 0.0
        T_bg = 1.0
        
        M_bg = HyQMOM.InitializeM4_35(rho_bg, u_bg, v_bg, w_bg, T_bg, 0.0, 0.0, T_bg, 0.0, T_bg)
        
        @test M_bg[1] ≈ rho_bg
        @test M_bg[2] ≈ 0.0 atol=1e-10  # M100
        @test M_bg[6] ≈ 0.0 atol=1e-10  # M010
        @test M_bg[16] ≈ 0.0 atol=1e-10 # M001
    end
    
    @testset "Multiple jets don't overlap" begin
        # Test that we can detect non-overlapping jets
        jet1_center = (0.3, 0.5, 0.5)
        jet2_center = (0.7, 0.5, 0.5)
        width = (0.2, 0.2, 0.2)
        
        jet1 = HyQMOM.CubicRegion(
            center=jet1_center,
            width=width,
            density=2.0,
            velocity=(10.0, 0.0, 0.0),
            temperature=1.0
        )
        
        jet2 = HyQMOM.CubicRegion(
            center=jet2_center,
            width=width,
            density=2.0,
            velocity=(-10.0, 0.0, 0.0),
            temperature=1.0
        )
        
        # A point in jet1 should not be in jet2
        point_in_jet1 = jet1_center
        @test HyQMOM.point_in_cube(point_in_jet1, jet1)
        @test !HyQMOM.point_in_cube(point_in_jet1, jet2)
        
        # A point in jet2 should not be in jet1
        point_in_jet2 = jet2_center
        @test HyQMOM.point_in_cube(point_in_jet2, jet2)
        @test !HyQMOM.point_in_cube(point_in_jet2, jet1)
    end
    
    @testset "Crossing jets geometry" begin
        # Test classic crossing jets setup
        offset = 0.12  # Jets separated by 2*offset
        jet_width = (0.2, 0.2, 0.2)
        
        # Horizontal jet (moving right)
        jet_x = HyQMOM.CubicRegion(
            center=(0.5 - offset, 0.5, 0.5),
            width=jet_width,
            density=2.0,
            velocity=(50.0, 0.0, 0.0),
            temperature=1.0
        )
        
        # Vertical jet (moving up)
        jet_y = HyQMOM.CubicRegion(
            center=(0.5, 0.5 - offset, 0.5),
            width=jet_width,
            density=2.0,
            velocity=(0.0, 50.0, 0.0),
            temperature=1.0
        )
        
        # Jets should not overlap at center
        center_point = (0.5, 0.5, 0.5)
        in_jet_x = HyQMOM.point_in_cube(center_point, jet_x)
        in_jet_y = HyQMOM.point_in_cube(center_point, jet_y)
        
        # For offset > width/2, neither should contain center
        if offset > jet_width[1]/2
            @test !in_jet_x
            @test !in_jet_y
        end
    end
    
    @testset "Temperature consistency" begin
        # Test that temperature is correctly initialized
        rho = 2.0
        u, v, w = 5.0, 3.0, -2.0
        T = 1.5
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # Extract temperature from second moments
        M200 = M[3]
        Txx = M200 / rho - u^2
        
        @test Txx ≈ T rtol=1e-10
    end
    
    @testset "Velocity consistency" begin
        # Test that velocity is correctly initialized
        rho = 2.0
        u, v, w = 10.0, -5.0, 3.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # Extract velocities from first moments
        M100 = M[2]
        M010 = M[6]
        M001 = M[16]
        
        u_extracted = M100 / rho
        v_extracted = M010 / rho
        w_extracted = M001 / rho
        
        @test u_extracted ≈ u rtol=1e-10
        @test v_extracted ≈ v rtol=1e-10
        @test w_extracted ≈ w rtol=1e-10
    end
    
    @testset "Symmetry in initial conditions" begin
        # Symmetric jets should produce symmetric moment fields
        rho = 2.0
        u = 50.0
        T = 1.0
        width = (0.2, 0.2, 0.2)
        
        # Symmetric setup: two jets with opposite velocities
        jet1 = HyQMOM.CubicRegion(
            center=(0.3, 0.5, 0.5),
            width=width,
            density=rho,
            velocity=(u, 0.0, 0.0),
            temperature=T
        )
        
        jet2 = HyQMOM.CubicRegion(
            center=(0.7, 0.5, 0.5),
            width=width,
            density=rho,
            velocity=(-u, 0.0, 0.0),
            temperature=T
        )
        
        # Both jets should have same density and temperature
        @test jet1.density == jet2.density
        @test jet1.temperature == jet2.temperature
        
        # Velocities should be opposite
        @test jet1.velocity[1] == -jet2.velocity[1]
    end
end

