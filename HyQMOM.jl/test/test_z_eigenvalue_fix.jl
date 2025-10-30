"""
Regression test for z-direction eigenvalue bug fix.

This test verifies that the z-eigenvalue moment ordering bug (commit 83288b6) 
stays fixed. The bug caused pathological eigenvalues (~3e5) due to incorrect 
moment ordering in the WV (z-y plane) Jacobian.

Key test: For Ma=70 crossing jets, eigenvalues should be ~50, not ~3e5.
"""

using Test
using HyQMOM

@testset "Z-eigenvalue bug regression test" begin
    
    @testset "WV Jacobian uses correct moment ordering" begin
        # Create moments that would trigger the bug
        # Use more realistic moments from Maxwell-Boltzmann distribution
        
        M = zeros(35)
        rho = 1.0
        u, v, w = 40.0, 40.0, -40.0
        T = 1.0
        
        M[1] = rho
        M[2] = rho * u
        M[6] = rho * v  
        M[16] = rho * w
        
        # Second moments: <v_i^2> = u_i^2 + T
        M[3] = rho * (u^2 + T)      # M200
        M[10] = rho * (v^2 + T)     # M020
        M[20] = rho * (w^2 + T)     # M002
        
        # Third moments: <v_i^3> = u_i^3 + 3*u_i*T
        M[4] = rho * (u^3 + 3*u*T)  # M300
        M[13] = rho * (v^3 + 3*v*T) # M030
        M[23] = rho * (w^3 + 3*w*T) # M003
        
        # Fourth moments: <v_i^4> = u_i^4 + 6*u_i^2*T + 3*T^2
        M[5] = rho * (u^4 + 6*u^2*T + 3*T^2)   # M400
        M[15] = rho * (v^4 + 6*v^2*T + 3*T^2)  # M040
        M[25] = rho * (w^4 + 6*w^2*T + 3*T^2)  # M004
        
        # Mixed moments: <v_i*v_j> = u_i*u_j (uncorrelated)
        M[7] = rho * u * v      # M110
        M[17] = rho * u * w     # M101
        M[26] = rho * v * w     # M011
        
        flag2D = 0
        Ma = 70.0
        
        # Compute z-direction eigenvalues
        v6zmin, v6zmax, _ = HyQMOM.eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
        
        # Before fix: v6zmax was ~3.3e5 (pathological)
        # After fix: v6zmax should be ~50 (physical)
        
        @test isfinite(v6zmin)
        @test isfinite(v6zmax)
        
        # Main regression test: eigenvalues must be physical
        @test abs(v6zmax) < 5000.0  # Way below pathological (3e5)
        @test abs(v6zmin) < 5000.0
        
        # Should be on order of velocity
        @test abs(v6zmax) < 5.0 * abs(Ma)  # Generous bound
    end
    
    @testset "All directions produce comparable eigenvalues" begin
        # For symmetric initial conditions, all three directions 
        # should produce similar eigenvalue magnitudes
        
        M = zeros(35)
        rho = 1.0
        u, v, w = 40.0, 40.0, -40.0
        T = 1.0
        
        M[1] = rho
        M[2], M[6], M[16] = rho*u, rho*v, rho*w
        
        # Second moments
        M[3] = rho * (u^2 + T)
        M[10] = rho * (v^2 + T)
        M[20] = rho * (w^2 + T)
        
        # Fourth moments (minimal)
        M[5] = rho * (u^4 + 6*u^2*T + 3*T^2)
        M[15] = rho * (v^4 + 6*v^2*T + 3*T^2)
        M[25] = rho * (w^4 + 6*w^2*T + 3*T^2)
        
        # Mixed moments
        M[7] = rho * u * v      # M110
        M[17] = rho * u * w     # M101
        M[26] = rho * v * w     # M011
        
        flag2D = 0
        Ma = 70.0
        
        # X-direction
        v6xmin, v6xmax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        
        # Y-direction
        v6ymin, v6ymax, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 2, flag2D, Ma)
        
        # Z-direction (this had the bug)
        v6zmin, v6zmax, _ = HyQMOM.eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
        
        # All finite
        @test isfinite(v6xmax) && isfinite(v6ymax) && isfinite(v6zmax)
        
        # All physical (not pathological)
        @test abs(v6xmax) < 5000.0
        @test abs(v6ymax) < 5000.0
        @test abs(v6zmax) < 5000.0  # Was 3.3e5 before fix!
        
        # All three directions comparable (within factor of 10)
        max_eig = max(abs(v6xmax), abs(v6ymax), abs(v6zmax))
        min_eig = min(abs(v6xmax), abs(v6ymax), abs(v6zmax))
        
        if min_eig > 1.0  # Only check ratio if eigenvalues are significant
            @test max_eig / min_eig < 10.0
        end
    end
    
    @testset "High Ma remains stable" begin
        # Test that fix works even at very high Mach numbers
        
        M = zeros(35)
        Ma = 150.0
        u = Ma / sqrt(3.0)
        T = 1.0
        
        M[1] = 1.0
        M[2], M[6], M[16] = u, u, -u
        
        M[3] = u^2 + T
        M[10] = u^2 + T
        M[20] = u^2 + T
        
        M[5] = (u^2 + T)^2
        M[15] = (u^2 + T)^2
        M[25] = (u^2 + T)^2
        
        flag2D = 0
        
        v6zmin, v6zmax, _ = HyQMOM.eigenvalues6z_hyperbolic_3D(M, flag2D, Ma)
        
        # Should not be pathological even at high Ma
        @test isfinite(v6zmax)
        @test abs(v6zmax) < 1e4  # Much less than 3e5
        
        # Should scale reasonably with Ma
        @test abs(v6zmax) < 50.0 * Ma  # Allow generous bound for high Ma
    end
end

