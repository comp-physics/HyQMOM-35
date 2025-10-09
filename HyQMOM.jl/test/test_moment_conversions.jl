const TOL = 1e-10

@testset "Moment Conversions" begin
    
    @testset "M2CS4_35 Gaussian" begin
        rho = 1.0
        u, v, w = 0.1, 0.2, 0.3
        T = 1.5
        
        M = InitializeM4_35(rho, u, v, w, T, 0.0, 0.0, T, 0.0, T)
        
        C4, S4 = M2CS4_35(M)
        
        @test C4[3] ~= T atol=TOL  # C200 should equal T
        @test C4[10] ~= T atol=TOL  # C020 should equal T
        @test C4[20] ~= T atol=TOL  # C002 should equal T
        
        @test S4[5] ~= 3.0 atol=TOL  # S400 should be 3 for Gaussian
        @test S4[15] ~= 3.0 atol=TOL  # S040 should be 3 for Gaussian
        @test S4[25] ~= 3.0 atol=TOL  # S004 should be 3 for Gaussian
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
end
