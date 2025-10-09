const TOL = 1e-10

@testset "Realizability" begin
    
    @testset "realizability_S2 basic" begin
        S110, S101, S011 = 0.5, 0.3, 0.4
        
        S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
        
        @test abs(S110r) <= 1.0 + TOL
        @test abs(S101r) <= 1.0 + TOL
        @test abs(S011r) <= 1.0 + TOL
        @test S2r >= 0.0
    end
    
    @testset "realizability_S111 basic" begin
        S110, S101, S011 = 0.5, 0.3, 0.4
        S210, S201, S120, S021, S102, S012 = 0.1, 0.2, 0.15, 0.25, 0.18, 0.22
        S111 = 0.3
        
        S111r = realizability(:S111, S110, S101, S011, S210, S201, S120, S021, S102, S012, S111)
        
        @test isfinite(S111r)
    end
    
    @testset "realizability_S210 basic" begin
        S110, S101, S011 = 0.5, 0.3, 0.4
        S300 = 0.1
        S210, S201 = 0.1, 0.2
        S400 = 3.0
        H200 = max(eps(), S400 - S300^2 - 1)
        beta = 1.0
        
        S210r, S201r = realizability(:S210, S110, S101, S011, S300, S210, S201, H200, beta)
        
        @test isfinite(S210r)
        @test isfinite(S201r)
    end
    
    @testset "realizability_S220 basic" begin
        S110 = 0.5
        S220 = 1.0
        A220 = 1.5
        
        S220r = realizability(:S220, S110, S220, A220)
        
        @test isfinite(S220r)
        @test abs(S220r) <= A220 + TOL
    end
    
    @testset "realizable_2D basic" begin
        S300, S400 = 0.1, 3.0
        S110, S210, S310 = 0.5, 0.1, 0.2
        S120, S220 = 0.15, 1.0
        S030, S130, S040 = 0.2, 0.25, 3.0
        
        S210r, S120r, S310r, S220r, S130r = realizability(Symbol("2D"), 
            S300, S400, S110, S210, S310, S120, S220, S030, S130, S040)
        
        @test all(isfinite.([S210r, S120r, S310r, S220r, S130r]))
    end
    
    @testset "Gaussian moments remain realizable" begin
        rho = 1.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.0
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        C4, S4 = M2CS4_35(M)
        
        # Extract 2nd-order moments
        S110, S101, S011 = S4[7], S4[17], S4[26]
        
        # Apply realizability
        S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
        
        # Should not change for Gaussian
        @test S110r ~= S110 atol=TOL
        @test S101r ~= S101 atol=TOL
        @test S011r ~= S011 atol=TOL
    end
    
    @testset "Realizability preserves finiteness" begin
        # Test with potentially problematic values
        S110, S101, S011 = 0.99, 0.98, 0.97
        
        S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
        
        @test all(isfinite.([S110r, S101r, S011r, S2r]))
        @test abs(S110r) <= 1.0 + TOL
        @test abs(S101r) <= 1.0 + TOL
        @test abs(S011r) <= 1.0 + TOL
    end
end
