const TOL = 1e-10

@testset "Realizability" begin
    
    @testset "realizability_S2 basic" begin
        S110, S101, S011 = 0.5, 0.3, 0.4
        
        S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
        
        @test abs(S110r) <= 1.0 + TOL
        @test abs(S101r) <= 1.0 + TOL
        @test abs(S011r) <= 1.0 + TOL
        @test S2r >= 0.0
        
        # Test with zero correlations
        S110r2, S101r2, S011r2, S2r2 = realizability(:S2, 0.0, 0.0, 0.0)
        @test S110r2 == 0.0
        @test S101r2 == 0.0
        @test S011r2 == 0.0
        @test all(isfinite.([S2r2]))
        
        # Test with extreme values (function handles clamping internally)
        S110r3, S101r3, S011r3, S2r3 = realizability(:S2, 1.5, 1.5, 1.5)
        @test all(isfinite.([S110r3, S101r3, S011r3, S2r3]))
        
        # Test with negative correlations
        S110r4, S101r4, S011r4, S2r4 = realizability(:S2, -0.8, -0.6, -0.7)
        @test all(isfinite.([S110r4, S101r4, S011r4, S2r4]))
        @test S2r4 >= 0.0
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
        
        # Test with zero correlations
        S210r2, S120r2, S310r2, S220r2, S130r2 = realizability(Symbol("2D"),
            0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 3.0)
        @test all(isfinite.([S210r2, S120r2, S310r2, S220r2, S130r2]))
        
        # Test realizable_2D directly (takes 10 parameters, returns 5: S21, S12, S31, S22, S13)
        result = HyQMOM.realizable_2D(
            S300, S400, S110, S210, S310, S120, S220,
            S030, S130, S040
        )
        
        @test length(result) == 5
        @test all(isfinite, result)
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
        @test S110r ≈ S110 atol=TOL
        @test S101r ≈ S101 atol=TOL
        @test S011r ≈ S011 atol=TOL
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

    @testset "realizable_3D direct call" begin
        # Test the core realizable_3D function directly with valid inputs
        # Function signature: realizable_3D(S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
        #                                   S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
        #                                   S211, S021, S121, S031, S012, S112, S013, S022)
        S300, S400 = 0.0, 3.0
        S110, S210, S310 = 0.0, 0.0, 0.0
        S120, S220 = 0.0, 1.5
        S030, S130, S040 = 0.0, 0.0, 3.0
        S101, S201, S301 = 0.0, 0.0, 0.0
        S102, S202 = 0.0, 0.0
        S003, S103, S004 = 0.0, 0.0, 3.0
        S011, S111, S211 = 0.0, 0.0, 0.0
        S021, S121, S031 = 0.0, 0.0, 0.0
        S012, S112, S013 = 0.0, 0.0, 0.0
        S022 = 1.5
        
        result = HyQMOM.realizable_3D(
            S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
            S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
            S211, S021, S121, S031, S012, S112, S013, S022
        )
        
        # Returns 29 values: all 28 standardized moments + flag220
        @test length(result) == 29
        @test all(isfinite, result)
        
        # Test with some correlations
        S110, S101, S011 = 0.3, 0.2, 0.1
        S111 = 0.05
        
        result2 = HyQMOM.realizable_3D(
            S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
            S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
            S211, S021, S121, S031, S012, S112, S013, S022
        )
        
        @test length(result2) == 29
        @test all(isfinite, result2)
    end
    
    @testset "realizability_3D basic" begin
        # Start from an isotropic Gaussian-like state: third-order and cross moments
        # are zero and kurtosis is ~3 in each direction.
        S300 = 0.0
        S400 = 3.0
        S030 = 0.0
        S040 = 3.0
        S003 = 0.0
        S004 = 3.0

        S110 = 0.0
        S101 = 0.0
        S011 = 0.0

        S210 = 0.0
        S310 = 0.0
        S120 = 0.0
        S220 = 0.0
        S130 = 0.0

        S201 = 0.0
        S301 = 0.0
        S102 = 0.0
        S103 = 0.0
        S202 = 0.0

        S111 = 0.0
        S211 = 0.0
        S021 = 0.0
        S121 = 0.0
        S031 = 0.0
        S012 = 0.0
        S112 = 0.0
        S013 = 0.0
        S022 = 0.0

        result = realizability(Symbol("3D"),
                               S300, S400, S110, S210, S310, S120, S220, S030, S130, S040,
                               S101, S201, S301, S102, S202, S003, S103, S004, S011, S111,
                               S211, S021, S121, S031, S012, S112, S013, S022)

        @test length(result) == 29

        vals = result[1:28]
        flag220 = result[29]
        @test all(isfinite.(vals))
        @test flag220 in (0, 1)
    end
    
    # Note: realizability_S211, S310, and S310_220 are tested indirectly through 
    # realizable_3D which calls them internally
    
    @testset "Edge cases - extreme correlations" begin
        # Test with high correlations (near 1)
        S110, S101, S011 = 0.99, 0.98, 0.97
        S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
        
        @test abs(S110r) <= 1.0 + TOL
        @test abs(S101r) <= 1.0 + TOL
        @test abs(S011r) <= 1.0 + TOL
        @test S2r >= 0.0
        
        # Test with negative correlations
        S110, S101, S011 = -0.5, -0.3, -0.4
        S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
        
        @test all(isfinite.([S110r, S101r, S011r, S2r]))
    end
    
    @testset "Realizability with varying temperatures" begin
        # Test that realizability works with different temperature values
        for T in [0.5, 1.0, 2.0, 5.0]
            rho = 1.0
            u, v, w = 0.5, 0.3, 0.0
            
            M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
            C4, S4 = M2CS4_35(M)
            
            # Extract and apply S2 realizability
            S110, S101, S011 = S4[7], S4[17], S4[26]
            S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
            
            @test all(isfinite.([S110r, S101r, S011r, S2r]))
        end
    end
    
    @testset "Realizability with varying velocities" begin
        rho = 1.0
        T = 1.0
        
        # Test with different velocity magnitudes
        for vel_mag in [0.0, 0.5, 1.0, 2.0, 5.0]
            u = vel_mag / sqrt(3.0)
            v = vel_mag / sqrt(3.0)
            w = vel_mag / sqrt(3.0)
            
            M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
            C4, S4 = M2CS4_35(M)
            
            S110, S101, S011 = S4[7], S4[17], S4[26]
            S110r, S101r, S011r, S2r = realizability(:S2, S110, S101, S011)
            
            @test all(isfinite.([S110r, S101r, S011r, S2r]))
        end
    end
    
    @testset "2D realizability with various inputs" begin
        # Test with different moment sets
        test_cases = [
            # (S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)
            (0.0, 3.0, 0.5, 0.1, 0.2, 0.15, 1.0, 0.0, 0.25, 3.0),
            (0.1, 3.5, 0.3, 0.05, 0.1, 0.08, 0.8, 0.05, 0.12, 3.2),
            (0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4.0),  # Uncorrelated
        ]
        
        for (S30, S40, S11, S21, S31, S12, S22, S03, S13, S04) in test_cases
            S21r, S12r, S31r, S22r, S13r = realizability(Symbol("2D"),
                S30, S40, S11, S21, S31, S12, S22, S03, S13, S04)
            
            @test all(isfinite.([S21r, S12r, S31r, S22r, S13r]))
        end
    end
    
    @testset "Realizability bounds checking" begin
        # Test S220 bounds
        S110 = 0.5
        A220 = 1.5
        
        # Test value that's too high (should be clamped to s220max = 1 + A220 = 2.5)
        S220_high = 3.0
        S220r_high = realizability(:S220, S110, S220_high, A220)
        s220max = 1.0 + A220
        @test S220r_high <= s220max + TOL
        @test S220r_high ≈ s220max atol=TOL
        
        # Test value that's too low (should be clamped to s220min = max(S110^2, 1-A220))
        S220_low = -1.0
        S220r_low = realizability(:S220, S110, S220_low, A220)
        s220min = max(S110^2, 1.0 - A220)
        @test S220r_low >= s220min - TOL
        @test S220r_low ≈ s220min atol=TOL
        
        # Test with S220 within bounds (should remain unchanged)
        S220_ok = 1.0
        S220r_ok = realizability(:S220, S110, S220_ok, A220)
        @test S220r_ok ≈ S220_ok atol=TOL
    end
end
