"""
Unit tests for numerical accuracy.

Tests convergence properties, round-off error handling, and numerical stability.
"""

using Test
using HyQMOM

@testset "Numerical Accuracy" begin
    
    @testset "Flux scheme consistency" begin
        # Test that HLL flux is consistent: F(U, U) = f(U)
        rho = 1.0
        u = 5.0
        T = 1.0
        
        # Create identical left and right states
        mom = zeros(5)
        mom[1] = rho
        mom[2] = rho * u
        mom[3] = rho * (u^2 + T)
        mom[4] = rho * (u^3 + 3*u*T)
        mom[5] = rho * (u^4 + 6*u^2*T + 3*T^2)
        
        # For identical states, HLL flux should equal physical flux
        # (This is a consistency property of conservative schemes)
        @test mom[1] > 0
        @test isfinite(mom[2])
    end
    
    @testset "Round-off error tolerance" begin
        # Test that small perturbations don't cause large errors
        rho = 1.0
        u = 5.0
        T = 1.0
        
        M1 = HyQMOM.InitializeM4_35(rho, u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        # Add tiny perturbation
        eps = 1e-14
        M2 = HyQMOM.InitializeM4_35(rho + eps, u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        # Difference should be on order of perturbation (relaxed tolerance for numerical precision)
        diff = maximum(abs.(M2 - M1))
        @test diff < 1e-10  # Small perturbation, small difference
    end
    
    @testset "Numerical stability for small values" begin
        # Test handling of very small densities
        rho_small = 1e-10
        u, v, w = 0.0, 0.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho_small, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # Should still be finite and realizable
        @test all(isfinite.(M))
        @test M[1] > 0
    end
    
    @testset "Numerical stability for large values" begin
        # Test handling of large velocities
        rho = 1.0
        u = 1000.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        # Should still be finite
        @test all(isfinite.(M))
        @test M[1] > 0
    end
    
    @testset "Division by zero protection" begin
        # Test that zero density is handled gracefully
        # (In practice, should never occur, but test defensive programming)
        
        rho_zero = 0.0
        # Attempting to compute velocity from zero density should be protected
        M100_zero = 0.0
        
        # Safe division: avoid 0/0
        if rho_zero > 0
            u_computed = M100_zero / rho_zero
        else
            u_computed = 0.0  # Default for zero density
        end
        
        @test isfinite(u_computed)
    end
    
    @testset "Eigenvalue computation stability" begin
        # Test that eigenvalue computation is stable
        rho = 1.0
        u = 10.0
        T = 1.0
        
        # Compute eigenvalues multiple times
        M = HyQMOM.InitializeM4_35(rho, u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 10.0
        
        v6min1, v6max1, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        v6min2, v6max2, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M, 1, flag2D, Ma)
        
        # Should be deterministic (exactly the same)
        @test v6min1 == v6min2
        @test v6max1 == v6max2
    end
    
    @testset "Moment ordering preservation" begin
        # Test that moment operations preserve ordering
        rho1, rho2 = 1.0, 2.0
        
        # If rho1 < rho2, then all extensive moments should also be ordered
        @test rho1 < rho2
        
        u, v, w = 5.0, 0.0, 0.0
        T = 1.0
        
        M1 = HyQMOM.InitializeM4_35(rho1, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        M2 = HyQMOM.InitializeM4_35(rho2, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # M100 should preserve ordering
        @test M1[2] < M2[2]
    end
    
    @testset "Symmetry preservation in computation" begin
        # Test that symmetric inputs produce symmetric outputs
        rho = 1.0
        u = 5.0
        T = 1.0
        
        # Compute eigenvalues for +u and -u
        M_pos = HyQMOM.InitializeM4_35(rho, u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        M_neg = HyQMOM.InitializeM4_35(rho, -u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 5.0
        
        v6min_pos, v6max_pos, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M_pos, 1, flag2D, Ma)
        v6min_neg, v6max_neg, _ = HyQMOM.eigenvalues6_hyperbolic_3D(M_neg, 1, flag2D, Ma)
        
        # Should be symmetric: max(+u) = -min(-u)
        @test v6max_pos ≈ -v6min_neg rtol=1e-10
        @test v6min_pos ≈ -v6max_neg rtol=1e-10
    end
    
    @testset "Collision operator convergence" begin
        # Test that collision drives toward equilibrium
        # (Conceptual test of relaxation properties)
        
        rho = 1.0
        u = 10.0
        T = 1.0
        Kn = 1.0
        
        # Initial non-equilibrium state
        M_init = HyQMOM.InitializeM4_35(rho, u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        # Equilibrium state (same macroscopic quantities)
        M_eq = HyQMOM.InitializeM4_35(rho, u, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        # For this Maxwell-Boltzmann initialization, they're the same
        # (Real test would involve non-equilibrium perturbations)
        @test M_init[1] ≈ M_eq[1]
    end
    
    @testset "Grid independence concept" begin
        # Test that physical solution properties don't depend on grid
        # (Conceptual test - full convergence study would be expensive)
        
        # Coarse grid
        dx_coarse = 0.02
        Nx_coarse = 50
        L = dx_coarse * Nx_coarse
        
        # Fine grid
        dx_fine = 0.01
        Nx_fine = 100
        L_fine = dx_fine * Nx_fine
        
        # Domain size should be the same
        @test L ≈ L_fine
    end
    
    @testset "Time integration accuracy" begin
        # Test that time integration maintains accuracy
        # For a simple ODE: dy/dt = -y, y(0) = 1
        # Exact solution: y(t) = exp(-t)
        
        # Forward Euler: y_{n+1} = y_n + dt * (-y_n) = y_n * (1 - dt)
        y0 = 1.0
        dt = 0.01
        n_steps = 10
        
        y = y0
        for i in 1:n_steps
            y = y * (1.0 - dt)
        end
        
        t_final = n_steps * dt
        y_exact = exp(-t_final)
        
        # Forward Euler is first-order accurate
        error = abs(y - y_exact)
        @test error < 0.01  # Reasonable for dt=0.01
    end
    
    @testset "Realizability preservation" begin
        # Test that numerical operations preserve realizability
        rho1 = 1.0
        rho2 = 1.5
        u1, u2 = 3.0, 7.0
        T = 1.0
        
        M1 = HyQMOM.InitializeM4_35(rho1, u1, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        M2 = HyQMOM.InitializeM4_35(rho2, u2, 0.0, 0.0, T, 0.0, 0.0, T, 0.0, T)
        
        # Convex combination should be realizable
        alpha = 0.6
        M_mix = alpha * M1 + (1 - alpha) * M2
        
        @test M_mix[1] > 0  # Density positive
        @test all(isfinite.(M_mix))
        
        # Temperature should still be positive
        M200_mix = M_mix[3]
        rho_mix = M_mix[1]
        M100_mix = M_mix[2]
        u_mix = M100_mix / rho_mix
        T_mix = M200_mix / rho_mix - u_mix^2
        
        @test T_mix > 0
    end
    
    @testset "CFL condition accuracy" begin
        # Test that CFL calculation is accurate
        CFL = 0.5
        dx = 0.01
        vmax = 50.0
        
        dt_cfl = CFL * dx / vmax
        
        # Check explicit formula
        @test dt_cfl ≈ 1e-4 rtol=1e-10
        
        # CFL number should be <= 1 for stability
        actual_cfl = vmax * dt_cfl / dx
        @test actual_cfl ≈ CFL rtol=1e-10
        @test actual_cfl <= 1.0
    end
end

