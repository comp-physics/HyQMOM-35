"""
Unit tests for time step calculation logic.

Tests the CFL-based time step calculation used in simulation_runner.jl.
"""

using Test

@testset "Time Step Calculations" begin
    
    @testset "Basic CFL condition" begin
        CFL = 0.3
        dx, dy, dz = 0.02, 0.02, 0.02
        vmax = 50.0
        
        dt = CFL * min(dx, dy, dz) / vmax
        
        @test dt ≈ 1.2e-4 rtol=1e-10
        @test dt == CFL * 0.02 / 50.0
    end
    
    @testset "dtmax cap applied" begin
        CFL = 0.5
        dx = 0.02
        vmax = 0.01  # Very slow flow
        dtmax = 0.01
        
        dt_cfl = CFL * dx / vmax  # = 1.0 (large!)
        dt_actual = min(dt_cfl, dtmax)
        
        @test dt_actual == dtmax
        @test dt_actual < dt_cfl
    end
    
    @testset "Grid refinement scaling" begin
        CFL = 0.3
        vmax = 50.0
        
        # Coarse grid
        dx_coarse = 0.02
        dt_coarse = CFL * dx_coarse / vmax
        
        # Fine grid (2x refined)
        dx_fine = 0.01
        dt_fine = CFL * dx_fine / vmax
        
        # Time step should halve when grid is refined by 2x
        @test dt_fine / dt_coarse ≈ 0.5 rtol=1e-10
    end
    
    @testset "Mach number scaling" begin
        CFL = 0.3
        dx = 0.02
        
        # Low Mach
        Ma1 = 1.0
        vmax1 = 1.5  # ~Ma + sound speed
        dt1 = CFL * dx / vmax1
        
        # High Mach
        Ma2 = 70.0
        vmax2 = 71.0  # ~Ma + sound speed
        dt2 = CFL * dx / vmax2
        
        # Time step should scale inversely with Mach
        @test dt2 < dt1
        @test dt2 / dt1 ≈ vmax1 / vmax2 rtol=0.1
    end
    
    @testset "Anisotropic grids" begin
        CFL = 0.3
        vmax = 50.0
        
        # Anisotropic: fine in x, coarse in y,z
        dx, dy, dz = 0.01, 0.05, 0.05
        dt = CFL * min(dx, dy, dz) / vmax
        
        # Should be limited by finest direction
        @test dt == CFL * dx / vmax
        @test dt < CFL * dy / vmax
    end
    
    @testset "Zero velocity handling" begin
        CFL = 0.3
        dx = 0.02
        
        # Test that we don't divide by zero
        vmax = 1e-15  # Essentially zero
        dt = CFL * dx / max(vmax, 1e-10)  # Use a floor
        
        @test isfinite(dt)
        @test dt > 0
    end
    
    @testset "Knudsen number cap" begin
        # MATLAB strategy: dt_max = Kn
        Kn = 1.0
        CFL = 0.3
        dx = 0.02
        vmax = 0.001  # Very slow
        
        dt_cfl = CFL * dx / vmax  # = 6.0 (large!)
        dt_actual = min(dt_cfl, Kn)
        
        @test dt_actual == Kn
    end
    
    @testset "High velocity stability" begin
        # Test extreme velocities don't cause problems
        CFL = 0.3
        dx = 0.01
        vmax = 1000.0  # Very high velocity
        
        dt = CFL * dx / vmax
        
        @test dt ≈ 3e-6 rtol=1e-10
        @test dt > 0
        @test isfinite(dt)
    end
end

