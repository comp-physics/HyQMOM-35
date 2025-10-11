const TOL = 1e-10

@testset "Initialization" begin
    
    @testset "InitializeM4_35 isotropic" begin
        rho = 1.5
        u, v, w = 0.1, 0.2, 0.3
        T = 2.0
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        @test M[1] ≈ rho atol=TOL
        @test M[2]/M[1] ≈ u atol=TOL
        @test M[6]/M[1] ≈ v atol=TOL
        @test M[16]/M[1] ≈ w atol=TOL
        
        C4, S4 = M2CS4_35(M)
        @test S4[5] ≈ 3.0 atol=TOL  # S400 kurtosis
        @test S4[15] ≈ 3.0 atol=TOL  # S040 kurtosis
        @test S4[25] ≈ 3.0 atol=TOL  # S004 kurtosis
    end
    
    @testset "InitializeM4_35 correlated" begin
        rho = 1.0
        u, v, w = 0.0, 0.0, 0.0
        C200, C020, C002 = 1.0, 1.5, 2.0
        C110, C101, C011 = 0.5, 0.3, 0.6
        
        M = InitializeM4_35(rho, u, v, w, C200, C110, C101, C020, C011, C002)
        
        C4, S4 = M2CS4_35(M)
        
        @test C4[3] ≈ C200 atol=TOL  # C200 variance
        @test C4[10] ≈ C020 atol=TOL  # C020 variance
        @test C4[20] ≈ C002 atol=TOL  # C002 variance
        
        S110_expected = C110 / sqrt(C200 * C020)
        @test S4[7] ≈ S110_expected atol=TOL  # S110 correlation
    end
    
    @testset "InitializeM4_35 mass conservation" begin
        rho = 3.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.5
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        @test M[1] ≈ rho atol=TOL
    end
    
    @testset "InitializeM4_35 third-order zero for Gaussian" begin
        rho = 1.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.0
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        C4, S4 = M2CS4_35(M)
        
        # Third-order standardized moments should be zero for Gaussian
        third_order_indices = [4, 8, 9, 11, 13, 18, 19, 21, 23, 24, 27, 28, 29, 30, 31, 32, 33, 34]
        
        for idx in third_order_indices
            @test abs(S4[idx]) < TOL
        end
    end
    
    @testset "InitializeM4_35 realizability" begin
        rho = 2.0
        u, v, w = 0.5, -0.3, 0.1
        T = 1.0
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # Basic realizability checks
        @test M[1] > 0  # M000 > 0
        @test all(isfinite.(M))
        
        # Check 2nd-order moments are positive definite
        C4, S4 = M2CS4_35(M)
        @test C4[3] > 0  # C200 > 0
        @test C4[10] > 0  # C020 > 0
        @test C4[20] > 0  # C002 > 0
    end
end
