"""
Unit tests for moment realizability bounds.

Tests physical constraints on moments beyond basic realizability checks.
"""

using Test
using HyQMOM

@testset "Realizability Bounds" begin
    
    @testset "Density positivity" begin
        # M000 (density) must always be positive
        rho = 1.5
        u, v, w = 10.0, 5.0, -3.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        @test M[1] > 0
        @test M[1] == rho
    end
    
    @testset "Temperature positivity from second moments" begin
        # Temperature inferred from second moments must be positive
        rho = 1.0
        u, v, w = 10.0, 5.0, -3.0
        T = 2.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # M200 = rho * (u^2 + Txx), so Txx = M200/rho - u^2
        M200 = M[3]
        Txx = M200 / rho - u^2
        @test Txx > 0
        @test Txx ≈ T rtol=1e-10
        
        # M020 = rho * (v^2 + Tyy)
        M020 = M[10]
        Tyy = M020 / rho - v^2
        @test Tyy > 0
        @test Tyy ≈ T rtol=1e-10
        
        # M002 = rho * (w^2 + Tzz)
        M002 = M[20]
        Tzz = M002 / rho - w^2
        @test Tzz > 0
        @test Tzz ≈ T rtol=1e-10
    end
    
    @testset "Cauchy-Schwarz for mixed moments" begin
        # M110^2 <= M200 * M020 (Cauchy-Schwarz inequality)
        rho = 1.0
        u, v, w = 10.0, 5.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        M200 = M[3]
        M020 = M[10]
        M110 = M[7]
        
        @test M110^2 <= M200 * M020
        
        # M101^2 <= M200 * M002
        M002 = M[20]
        M101 = M[17]
        @test M101^2 <= M200 * M002
        
        # M011^2 <= M020 * M002
        M011 = M[26]
        @test M011^2 <= M020 * M002
    end
    
    @testset "Higher moment bounds" begin
        # Third moments have bounds related to second moments
        rho = 1.0
        u, v, w = 5.0, 0.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        # All moments should be finite
        for i in 1:35
            @test isfinite(M[i])
        end
        
        # Fourth moments should be positive
        M400 = M[5]
        M040 = M[15]
        M004 = M[25]
        
        @test M400 > 0
        @test M040 > 0
        @test M004 > 0
    end
    
    @testset "Moment ordering" begin
        # For a distribution with u > 0, higher x-moments should increase
        rho = 1.0
        u = 10.0
        v, w = 0.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        M000 = M[1]   # rho
        M100 = M[2]   # rho * u
        M200 = M[3]   # rho * (u^2 + T)
        M300 = M[4]   # rho * (u^3 + 3*u*T)
        M400 = M[5]   # rho * (u^4 + 6*u^2*T + 3*T^2)
        
        # For u > 1, successive moments should generally increase
        @test M100 / M000 > 1.0  # Mean velocity > 0
        @test M200 / M100 > 1.0  # Second moment larger
    end
    
    @testset "Zero velocity realizability" begin
        # At rest (u=v=w=0), moments simplify
        rho = 1.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.0
        
        M = HyQMOM.InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        @test M[1] ≈ rho
        @test M[2] ≈ 0.0 atol=1e-10  # M100
        @test M[6] ≈ 0.0 atol=1e-10  # M010
        @test M[16] ≈ 0.0 atol=1e-10 # M001
        
        # Even moments are non-zero (temperature)
        @test M[3] ≈ rho * T  # M200
        @test M[10] ≈ rho * T # M020
        @test M[20] ≈ rho * T # M002
    end
    
    @testset "Realizability after closure" begin
        # Test that closure_and_eigenvalues maintains realizability
        rho = 1.0
        u = 5.0
        T = 1.0
        
        # 1D moments
        mom = zeros(5)
        mom[1] = rho
        mom[2] = rho * u
        mom[3] = rho * (u^2 + T)
        mom[4] = rho * (u^3 + 3*u*T)
        mom[5] = rho * (u^4 + 6*u^2*T + 3*T^2)
        
        Mp, vpmin, vpmax = HyQMOM.closure_and_eigenvalues(mom)
        
        @test Mp >= 0
        @test isfinite(Mp)
        @test isfinite(vpmin)
        @test isfinite(vpmax)
    end
end

