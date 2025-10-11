function tests = test_collision
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-10;
end

function test_collision35_mass_conservation(testCase)
    rho = 1.5;
    u = 0.1; v = 0.2; w = 0.3;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    dt = 0.01;
    Kn = 1.0;
    
    M_out = collision35(M, dt, Kn);
    
    verifyEqual(testCase, M_out(1), M(1), 'AbsTol', testCase.TestData.tol, ...
                'Mass should be conserved');
    
    fprintf('PASS: collision35 mass conservation test\n');
end

function test_collision35_momentum_conservation(testCase)
    rho = 1.0;
    u = 0.5; v = -0.3; w = 0.2;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    dt = 0.01;
    Kn = 1.0;
    
    M_out = collision35(M, dt, Kn);
    
    verifyEqual(testCase, M_out(2), M(2), 'AbsTol', testCase.TestData.tol, ...
                'X-momentum should be conserved');
    verifyEqual(testCase, M_out(6), M(6), 'AbsTol', testCase.TestData.tol, ...
                'Y-momentum should be conserved');
    verifyEqual(testCase, M_out(16), M(16), 'AbsTol', testCase.TestData.tol, ...
                'Z-momentum should be conserved');
    
    fprintf('PASS: collision35 momentum conservation test\n');
end

function test_collision35_relaxation_to_maxwellian(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    C200 = 1.5; C020 = 1.2; C002 = 1.8;
    C110 = 0.3; C101 = 0.2; C011 = 0.25;
    
    M = InitializeM4_35(rho, u, v, w, C200, C110, C101, C020, C011, C002);
    
    dt = 10.0;
    Kn = 0.01;
    
    M_out = collision35(M, dt, Kn);
    
    [~, S_out] = M2CS4_35(M_out);
    
    S110_out = S_out(7);
    S101_out = S_out(17);
    S011_out = S_out(26);
    
    verifyLessThan(testCase, abs(S110_out), abs(C110 / sqrt(C200 * C020)), ...
                   'S110 should relax toward zero');
    
    fprintf('PASS: collision35 relaxation to Maxwellian test\n');
end

function test_collision35_various_knudsen(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    dt = 0.01;
    Kn_values = [0.01, 0.1, 1.0, 10.0];
    
    for i = 1:length(Kn_values)
        Kn = Kn_values(i);
        M_out = collision35(M, dt, Kn);
        
        verifyEqual(testCase, M_out(1), M(1), 'AbsTol', testCase.TestData.tol, ...
                    sprintf('Mass conserved for Kn=%.2f', Kn));
        verifyTrue(testCase, all(isfinite(M_out)), ...
                   sprintf('All moments finite for Kn=%.2f', Kn));
    end
    
    fprintf('PASS: collision35 various Knudsen numbers test\n');
end

function test_collision35_temperature_relaxation(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    C200 = 2.0; C020 = 1.0; C002 = 1.5;
    
    M = InitializeM4_35(rho, u, v, w, C200, 0, 0, C020, 0, C002);
    
    dt = 1.0;
    Kn = 0.1;
    
    M_out = collision35(M, dt, Kn);
    
    [C_out, ~] = M2CS4_35(M_out);
    C200_out = C_out(3);
    C020_out = C_out(10);
    C002_out = C_out(20);
    
    T_initial = (C200 + C020 + C002) / 3;
    T_final = (C200_out + C020_out + C002_out) / 3;
    
    verifyEqual(testCase, T_final, T_initial, 'AbsTol', testCase.TestData.tol, ...
                'Average temperature should be conserved');
    
    variance_initial = var([C200, C020, C002]);
    variance_final = var([C200_out, C020_out, C002_out]);
    
    verifyLessThan(testCase, variance_final, variance_initial, ...
                   'Temperature components should become more isotropic');
    
    fprintf('PASS: collision35 temperature relaxation test\n');
end

function test_collision35_realizability(testCase)
    rho = 1.5;
    u = 0.2; v = -0.1; w = 0.3;
    T = 1.2;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    dt = 0.01;
    Kn = 1.0;
    
    M_out = collision35(M, dt, Kn);
    
    [realizable, violations] = verify_realizability(M_out, testCase.TestData.tol);
    
    verifyTrue(testCase, realizable, ...
               sprintf('Output moments should be realizable. Violations: %d', violations.count));
    
    fprintf('PASS: collision35 realizability test\n');
end
