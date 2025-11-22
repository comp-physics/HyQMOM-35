const TOL = 1e-10

@testset "Closures and Eigenvalues" begin
    
    @testset "closure_and_eigenvalues basic" begin
        # Create a simple moment vector
        mom = zeros(5)
        mom[1] = 1.0  # M0
        mom[2] = 0.1  # M1
        mom[3] = 1.2  # M2
        mom[4] = 0.15  # M3
        mom[5] = 3.5  # M4
        
        Mp, vpmin, vpmax = closure_and_eigenvalues(mom)
        
        @test isfinite(Mp)
        @test isfinite(vpmin)
        @test isfinite(vpmax)
        @test vpmin <= vpmax
    end
    
    @testset "Flux_closure35_and_realizable_3D produces output" begin
        rho = 1.0
        u, v, w = 0.1, 0.0, 0.0
        T = 1.0
        
        M4 = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        flag2D = 0
        Ma = 2.0
        
        # This should complete without errors
        Fx, Fy, Fz, M4r = Flux_closure35_and_realizable_3D(M4, flag2D, Ma)
        
        @test length(Fx) == 35
        @test length(Fy) == 35
        @test length(Fz) == 35
        @test length(M4r) == 35
        @test all(isfinite.(Fx))
        @test all(isfinite.(Fy))
        @test all(isfinite.(Fz))
        @test all(isfinite.(M4r))
    end
    
    @testset "Flux_closure35 preserves mass" begin
        rho = 2.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.5
        
        M4 = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        Fx, Fy, Fz, M4r = Flux_closure35_and_realizable_3D(M4, 0, 2.0)
        
        # M000 should be preserved
        @test M4r[1] ≈ M4[1] atol=TOL
    end
    
    @testset "eigenvalues6_hyperbolic_3D basic" begin
        rho = 1.0
        u, v, w = 0.1, 0.2, 0.3
        T = 1.0
        
        M4 = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # eigenvalues6_hyperbolic_3D(M::AbstractVector, axis::Int, flag2D::Int, Ma::Real)
        axis = 1  # X-axis
        flag2D = 0  # 3D
        Ma = 2.0
        
        vpmin, vpmax = eigenvalues6_hyperbolic_3D(M4, axis, flag2D, Ma)
        
        @test isfinite(vpmin)
        @test isfinite(vpmax)
        @test vpmin <= vpmax
    end
end
