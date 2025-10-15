const TOL = 1e-10

@testset "Numerical Schemes" begin
    
    @testset "flux_HLL basic" begin
        # flux_HLL(Wstar, W, l1, l2, F, N)
        # Simple test with grid of states
        N = 5
        Nmom = 3
        Wstar = ones(N, Nmom)
        W = ones(N, Nmom) * 1.1
        l1 = -ones(N)
        l2 = ones(N)
        F = zeros(N, Nmom)
        
        result = flux_HLL(Wstar, W, l1, l2, F, N)
        
        @test size(result) == (N-1, Nmom)
        @test all(isfinite.(result))
    end
    
    @testset "pas_HLL basic" begin
        Nx = 10
        Ny = 10
        Nmom = 5
        
        M = ones(Np, Nmom)
        F = zeros(Np, Nmom)
        
        # Simple flux
        for i in 1:Np
            F[i, :] = 0.1 * M[i, :]
        end
        
        dt = 0.01
        dx = 0.1
        vpmin = -ones(Np)
        vpmax = ones(Np)
        
        Mp = pas_HLL(M, F, dt, dx, vpmin, vpmax)
        
        @test size(Mp) == size(M)
        @test all(isfinite.(Mp))
    end
    
    @testset "pas_HLL with boundary conditions" begin
        Nx = 10
        Ny = 10
        Nmom = 5
        
        M = ones(Np, Nmom)
        F = zeros(Np, Nmom)
        
        dt = 0.01
        dx = 0.1
        vpmin = -ones(Np)
        vpmax = ones(Np)
        
        # Test with different BC combinations
        Mp1 = pas_HLL(M, F, dt, dx, vpmin, vpmax; apply_bc_left=true, apply_bc_right=true)
        Mp2 = pas_HLL(M, F, dt, dx, vpmin, vpmax; apply_bc_left=false, apply_bc_right=true)
        Mp3 = pas_HLL(M, F, dt, dx, vpmin, vpmax; apply_bc_left=true, apply_bc_right=false)
        
        @test all(isfinite.(Mp1))
        @test all(isfinite.(Mp2))
        @test all(isfinite.(Mp3))
    end
    
    @testset "collision35 basic" begin
        rho = 1.0
        u, v, w = 0.1, 0.2, 0.3
        T = 1.0
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        Kn = 0.01
        dt = 0.001
        
        Mnew = collision35(M, Kn, dt)
        
        @test length(Mnew) == 35
        @test all(isfinite.(Mnew))
        
        # Mass should be conserved
        @test Mnew[1] ≈ M[1] atol=TOL
    end
    
    @testset "collision35 relaxes to equilibrium" begin
        rho = 1.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.0
        
        # Start with non-equilibrium moments
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        M[11] += 0.1  # Perturb M300
        
        Kn = 0.01
        dt = 0.1  # Large time step
        
        Mnew = collision35(M, Kn, dt)
        
        # Should move toward equilibrium
        @test abs(Mnew[11] - M[11]) > 0
    end
end
