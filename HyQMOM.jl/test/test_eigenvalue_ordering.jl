"""
Unit tests for eigenvalue ordering and properties.

Tests that eigenvalues maintain correct ordering, symmetry, and physical bounds.
"""

using Test
using HyQMOM

@testset "Eigenvalue Ordering and Properties" begin
    
    @testset "Eigenvalues sorted correctly" begin
        # Eigenvalues should satisfy λ_min <= λ_max
        rho = 1.0
        u, v, w = 5.0, 3.0, -2.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 6.0
        
        # X-direction eigenvalues
        v6xmin, v6xmax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        @test v6xmin <= v6xmax
        
        # Y-direction eigenvalues
        v6ymin, v6ymax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 2, flag2D, Ma)
        @test v6ymin <= v6ymax
        
        # Z-direction eigenvalues
        v6zmin, v6zmax, _ = HyQMOM.eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
        @test v6zmin <= v6zmax
    end
    
    @testset "Eigenvalue symmetry for symmetric flows" begin
        # For symmetric velocities, all directions should have similar eigenvalues
        rho = 1.0
        u, v, w = 50.0, 50.0, 50.0  # Symmetric velocities
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = sqrt(3) * 50.0  # sqrt(u^2 + v^2 + w^2)
        
        v6xmin, v6xmax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        v6ymin, v6ymax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 2, flag2D, Ma)
        v6zmin, v6zmax, _ = HyQMOM.eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
        
        # Max eigenvalues should be comparable (within 10%)
        @test isapprox(v6xmax, v6ymax, rtol=0.1)
        @test isapprox(v6xmax, v6zmax, rtol=0.1)
        @test isapprox(v6ymax, v6zmax, rtol=0.1)
        
        # Min eigenvalues should be comparable
        @test isapprox(v6xmin, v6ymin, rtol=0.1)
        @test isapprox(v6xmin, v6zmin, rtol=0.1)
        @test isapprox(v6ymin, v6zmin, rtol=0.1)
    end
    
    @testset "Eigenvalues bounded by physical velocity" begin
        # Eigenvalues should be on order of flow velocity + sound speed
        rho = 1.0
        u = 10.0
        v, w = 0.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 10.0
        
        v6xmin, v6xmax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        
        # Sound speed ~ sqrt(T) ~ 1
        c_sound = sqrt(T)
        
        # Max eigenvalue should be roughly u + c_sound
        @test v6xmax > u
        @test v6xmax < 5.0 * (abs(u) + c_sound)  # Generous bound
        
        # Min eigenvalue should be roughly u - c_sound
        @test v6xmin < u
        @test v6xmin > -5.0 * (abs(u) + c_sound)
    end
    
    @testset "High Mach number stability" begin
        # Test that high Mach numbers don't cause pathological eigenvalues
        rho = 1.0
        u = 100.0
        v, w = 0.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 100.0
        
        v6xmin, v6xmax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        
        # Should be finite and not pathological
        @test isfinite(v6xmin)
        @test isfinite(v6xmax)
        @test abs(v6xmax) < 10000.0  # Not pathological (should be ~ 100)
        @test abs(v6xmin) < 10000.0
    end
    
    @testset "Eigenvalues for zero velocity" begin
        # At rest, eigenvalues should be small (sound speed only)
        rho = 1.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 0.0
        
        v6xmin, v6xmax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        
        # Should be symmetric around zero
        @test abs(v6xmin + v6xmax) < 1.0  # Sum should be near zero
        @test abs(v6xmax) < 5.0  # Should be on order of sound speed
    end
    
    @testset "Eigenvalue sign consistency" begin
        # Positive velocity should give positive max eigenvalue
        rho = 1.0
        u = 20.0
        v, w = 0.0, 0.0
        T = 1.0
        
        M_pos = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 20.0
        
        v6xmin_pos, v6xmax_pos, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M_pos, 1, flag2D, Ma)
        
        @test v6xmax_pos > 0
        
        # Negative velocity should give negative min eigenvalue
        u_neg = -20.0
        M_neg = HyQMOM.InitializeM4_35(rho, u_neg, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        v6xmin_neg, v6xmax_neg, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M_neg, 1, flag2D, Ma)
        
        @test v6xmin_neg < 0
    end
    
    @testset "Z-eigenvalue consistency with X and Y" begin
        # Z-eigenvalues should follow same patterns as X and Y
        rho = 1.0
        u, v, w = 10.0, 5.0, 15.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = sqrt(u^2 + v^2 + w^2)
        
        v6xmin, v6xmax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        v6ymin, v6ymax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 2, flag2D, Ma)
        v6zmin, v6zmax, _ = HyQMOM.eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
        
        # All should be finite
        @test isfinite(v6xmax) && isfinite(v6ymax) && isfinite(v6zmax)
        @test isfinite(v6xmin) && isfinite(v6ymin) && isfinite(v6zmin)
        
        # All should be on same order of magnitude
        max_eig = maximum([abs(v6xmax), abs(v6ymax), abs(v6zmax)])
        min_eig = minimum([abs(v6xmax), abs(v6ymax), abs(v6zmax)])
        
        # Should not differ by more than a factor of 10 for this case
        @test max_eig / min_eig < 10.0
    end
    
    @testset "1D closure eigenvalues" begin
        # Test 1D closure function eigenvalues
        rho = 1.0
        u = 5.0
        T = 1.0
        
        mom = zeros(5)
        mom[1] = rho
        mom[2] = rho * u
        mom[3] = rho * (u^2 + T)
        mom[4] = rho * (u^3 + 3*u*T)
        mom[5] = rho * (u^4 + 6*u^2*T + 3*T^2)
        
        Mp, vpmin, vpmax = HyQMOM.closure_and_eigenvalues(mom)
        
        @test vpmin <= vpmax
        @test isfinite(vpmin)
        @test isfinite(vpmax)
        
        # Should be bounded by velocity
        @test vpmax > u
        @test vpmin < u
    end
    
    @testset "Eigenvalue regression: no pathological values" begin
        # Regression test for z-eigenvalue bug
        # This specifically tests the bug fixed in eigenvalues6z_hyperbolic_3D
        rho = 1.0
        u, v, w = 40.0, 40.0, -40.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 70.0
        
        v6zmin, v6zmax, _ = HyQMOM.eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
        
        # Before fix: v6zmax was ~3.3e5 (pathological)
        # After fix: v6zmax should be ~50-100 (physical)
        @test abs(v6zmax) < 5000.0  # Well below pathological
        @test abs(v6zmin) < 5000.0
        
        # Should be on order of velocity + sound speed
        @test abs(v6zmax) < 5.0 * Ma
    end
end

