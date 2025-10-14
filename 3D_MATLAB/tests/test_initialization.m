function tests = test_initialization
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-12;
end

function test_InitializeM4_35_isotropic(testCase)
    rho = 1.5;
    u = 0.1; v = 0.2; w = 0.3;
    T = 2.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    verifyEqual(testCase, M(1), rho, 'AbsTol', testCase.TestData.tol, 'M000 should equal rho');
    verifyEqual(testCase, M(2)/M(1), u, 'AbsTol', testCase.TestData.tol, 'Mean velocity u');
    verifyEqual(testCase, M(6)/M(1), v, 'AbsTol', testCase.TestData.tol, 'Mean velocity v');
    verifyEqual(testCase, M(16)/M(1), w, 'AbsTol', testCase.TestData.tol, 'Mean velocity w');
    
    [~, S4] = M2CS4_35(M);
    verifyEqual(testCase, S4(5), 3.0, 'AbsTol', testCase.TestData.tol, 'S400 kurtosis');
    verifyEqual(testCase, S4(15), 3.0, 'AbsTol', testCase.TestData.tol, 'S040 kurtosis');
    verifyEqual(testCase, S4(25), 3.0, 'AbsTol', testCase.TestData.tol, 'S004 kurtosis');
    
    fprintf('PASS: InitializeM4_35 isotropic test\n');
end

function test_InitializeM4_35_correlated(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    C200 = 1.0; C020 = 1.5; C002 = 2.0;
    C110 = 0.5; C101 = 0.3; C011 = 0.6;
    
    M = InitializeM4_35(rho, u, v, w, C200, C110, C101, C020, C011, C002);
    
    [C4, S4] = M2CS4_35(M);
    
    verifyEqual(testCase, C4(3), C200, 'AbsTol', testCase.TestData.tol, 'C200 variance');
    verifyEqual(testCase, C4(10), C020, 'AbsTol', testCase.TestData.tol, 'C020 variance');
    verifyEqual(testCase, C4(20), C002, 'AbsTol', testCase.TestData.tol, 'C002 variance');
    
    S110_expected = C110 / sqrt(C200 * C020);
    verifyEqual(testCase, S4(7), S110_expected, 'AbsTol', testCase.TestData.tol, 'S110 correlation');
    
    fprintf('PASS: InitializeM4_35 correlated test\n');
end

function test_InitializeM4_35_realizability(testCase)
    rho = 2.0;
    u = 0.5; v = -0.3; w = 0.1;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [realizable, violations] = verify_realizability(M, testCase.TestData.tol);
    
    verifyTrue(testCase, realizable, ...
               sprintf('Initialized moments should be realizable. Violations: %d', violations.count));
    
    fprintf('PASS: InitializeM4_35 realizability test\n');
end

function test_InitializeM4_35_mass_conservation(testCase)
    rho = 3.0;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.5;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    verifyEqual(testCase, M(1), rho, 'AbsTol', testCase.TestData.tol, ...
                'Mass (M000) should equal input density');
    
    fprintf('PASS: InitializeM4_35 mass conservation test\n');
end

function test_InitializeM4_35_third_order_zero(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    [~, S4] = M2CS4_35(M);
    
    third_order_indices = [4, 8, 9, 11, 13, 18, 19, 21, 23, 24, 27, 28, 29, 30, 31, 32, 33, 34];
    
    for idx = third_order_indices
        verifyEqual(testCase, S4(idx), 0.0, 'AbsTol', testCase.TestData.tol, ...
                    sprintf('Third-order moment S4(%d) should be zero for Gaussian', idx));
    end
    
    fprintf('PASS: InitializeM4_35 third-order moments zero test\n');
end
