const TOL = 1e-10

@testset "Moment Conversions" begin
    
    @testset "M2CS4_35 Gaussian" begin
        rho = 1.0
        u, v, w = 0.1, 0.2, 0.3
        T = 1.5
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        C4, S4 = M2CS4_35(M)
        
        @test C4[3] ≈ T atol=TOL  # C200 should equal T
        @test C4[10] ≈ T atol=TOL  # C020 should equal T
        @test C4[20] ≈ T atol=TOL  # C002 should equal T
        
        @test S4[5] ≈ 3.0 atol=TOL  # S400 should be 3 for Gaussian
        @test S4[15] ≈ 3.0 atol=TOL  # S040 should be 3 for Gaussian
        @test S4[25] ≈ 3.0 atol=TOL  # S004 should be 3 for Gaussian
    end
    
    @testset "M4toC4_3D produces correct output" begin
        rho = 2.0
        u, v, w = 0.5, -0.3, 0.1
        C200, C020, C002 = 1.2, 1.5, 1.8
        C110, C101, C011 = 0.3, -0.2, 0.4
        
        M_orig = InitializeM4_35(rho, u, v, w, C200, C110, C101, C020, C011, C002)
        
        M000, M100, M010, M001, M200, M110, M101, M020, M011, M002,
        M300, M210, M201, M120, M111, M102, M030, M021, M012, M003,
        M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004 =
            M4_to_vars(M_orig)
        
        C_array = M4toC4_3D(M000, M100, M010, M001, M200, M110, M101, M020, M011, M002,
                            M300, M210, M201, M120, M111, M102, M030, M021, M012, M003,
                            M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004)
        
        @test size(C_array) == (5, 5, 5)
        @test all(isfinite.(C_array))
    end
    
    @testset "Standardized moments bounds" begin
        rho = 1.0
        u, v, w = 0.0, 0.0, 0.0
        T = 1.0
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        C4, S4 = M2CS4_35(M)
        
        S110 = S4[7]
        S101 = S4[17]
        S011 = S4[26]
        
        @test abs(S110) <= 1.0 + TOL
        @test abs(S101) <= 1.0 + TOL
        @test abs(S011) <= 1.0 + TOL
    end
    
    @testset "moment_idx" begin
        @test moment_idx("M110") == 7
        @test moment_idx("M200") == 3
        @test moment_idx("M000") == 1
        @test moment_idx("M100") == 2
        @test moment_idx("M010") == 6
        @test moment_idx("M001") == 16
    end
    
    @testset "M4_to_vars vector" begin
        M = rand(35)
        M[1] = abs(M[1]) + 1.0
        
        vars = M4_to_vars(M)
        
        @test length(vars) == 35
        @test vars[1] == M[1]
        @test vars[7] == M[7]
        @test vars[35] == M[35]
    end
    
    @testset "M4_to_vars 3D array" begin
        M_array = rand(5, 5, 5)
        
        vars = M4_to_vars(M_array)
        
        @test length(vars) == 35
        @test vars[1] == M_array[1,1,1]  # M000
    end

    @testset "moment_idx order 4 and names consistency" begin
        idxs = moment_idx(4)
        names = HyQMOM.moment_names(4)

        @test length(idxs) == 35
        @test length(names) == 35

        for (pos, name) in enumerate(names)
            i = parse(Int, name[2:2]) + 1
            j = parse(Int, name[3:3]) + 1
            k = parse(Int, name[4:4]) + 1
            linear_idx = i + 5*(j-1) + 25*(k-1)
            @test idxs[pos] == linear_idx
            @test moment_idx(name) == pos
        end

        @test_throws ErrorException moment_idx(5)
    end

    @testset "M5_to_vars vector" begin
        M5 = collect(1.0:56.0)
        vars = M5_to_vars(M5)

        @test length(vars) == 56
        @test vars[1] == M5[1]
        @test vars[20] == M5[20]
        @test vars[56] == M5[56]
    end

    @testset "M5_to_vars 3D array" begin
        M5 = zeros(6, 6, 6)
        M5[1,1,1] = 1.0  # M000
        M5[2,1,1] = 2.0  # M100
        M5[1,2,1] = 3.0  # M010
        M5[1,1,2] = 4.0  # M001

        vars = M5_to_vars(M5)

        @test length(vars) == 56
        @test vars[1] == 1.0   # M000
        @test vars[2] == 2.0   # M100
        @test vars[3] == 3.0   # M010
        @test vars[4] == 4.0   # M001
    end

    @testset "S_to_C_batch basic scaling" begin
        S110 = 0.5
        S101 = 0.25
        S011 = -0.3

        # Set a few standardized moments nonzero, others zero
        S300 = 1.0
        S210 = 0.0
        S201 = 0.0
        S120 = 0.0
        S111 = 0.0
        S102 = 0.0
        S030 = 2.0
        S021 = 0.0
        S012 = 0.0
        S003 = 3.0

        S400 = 0.0
        S310 = 0.0
        S301 = 0.0
        S220 = 0.0
        S211 = 0.0
        S202 = 0.0
        S130 = 0.0
        S121 = 0.0
        S112 = 0.0
        S103 = 0.0
        S040 = 0.0
        S031 = 0.0
        S022 = 0.0
        S013 = 0.0
        S004 = 0.0

        S500 = 0.0
        S410 = 0.0
        S401 = 0.0
        S320 = 0.0
        S311 = 0.0
        S302 = 0.0
        S230 = 0.0
        S221 = 0.0
        S212 = 0.0
        S203 = 0.0
        S140 = 0.0
        S131 = 0.0
        S122 = 0.0
        S113 = 0.0
        S104 = 0.0
        S050 = 0.0
        S041 = 0.0
        S032 = 0.0
        S023 = 0.0
        S014 = 0.0
        S005 = 0.0

        sC200, sC020, sC002 = 2.0, 3.0, 5.0

        C = HyQMOM.S_to_C_batch(S110, S101, S011, S300, S210, S201, S120, S111, S102, S030, S021, S012, S003,
                                 S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004,
                                 S500, S410, S401, S320, S311, S302, S230, S221, S212, S203, S140, S131, S122, S113, S104,
                                 S050, S041, S032, S023, S014, S005,
                                 sC200, sC020, sC002)

        @test length(C) ≥ 13

        C110 = C[1]
        C101 = C[2]
        C011 = C[3]
        C300 = C[4]
        C030 = C[10]
        C003 = C[13]

        @test C110 ≈ S110 * sC200 * sC020 atol=TOL
        @test C101 ≈ S101 * sC200 * sC002 atol=TOL
        @test C011 ≈ S011 * sC020 * sC002 atol=TOL
        @test C300 ≈ S300 * sC200^3 atol=TOL
        @test C030 ≈ S030 * sC020^3 atol=TOL
        @test C003 ≈ S003 * sC002^3 atol=TOL
    end
end
