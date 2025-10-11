function tests = test_moment_conversions
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-12;
    testCase.TestData.script_dir = script_dir;
end

function test_M2CS4_35_gaussian(testCase)
    rho = 1.0;
    u = 0.1; v = 0.2; w = 0.3;
    T = 1.5;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [C4, S4] = M2CS4_35(M);
    
    verifyEqual(testCase, C4(3), T, 'AbsTol', testCase.TestData.tol, 'C200 should equal T');
    verifyEqual(testCase, C4(10), T, 'AbsTol', testCase.TestData.tol, 'C020 should equal T');
    verifyEqual(testCase, C4(20), T, 'AbsTol', testCase.TestData.tol, 'C002 should equal T');
    
    verifyEqual(testCase, S4(5), 3.0, 'AbsTol', testCase.TestData.tol, 'S400 should be 3 for Gaussian');
    verifyEqual(testCase, S4(15), 3.0, 'AbsTol', testCase.TestData.tol, 'S040 should be 3 for Gaussian');
    verifyEqual(testCase, S4(25), 3.0, 'AbsTol', testCase.TestData.tol, 'S004 should be 3 for Gaussian');
    
    fprintf('PASS: M2CS4_35 Gaussian test\n');
end

function test_moment_struct_roundtrip(testCase)
    M_vec = rand(35, 1);
    M_vec(1) = abs(M_vec(1)) + 1.0;
    
    M_struct = moment_struct('from_vector', M_vec);
    M_vec2 = moment_struct('to_vector', M_struct);
    
    verifyEqual(testCase, M_vec, M_vec2, 'AbsTol', testCase.TestData.tol, ...
                'Vector-struct-vector roundtrip should preserve values');
    
    fprintf('PASS: moment_struct roundtrip test\n');
end

function test_M4toC4_3D_produces_output(testCase)
    rho = 2.0;
    u = 0.5; v = -0.3; w = 0.1;
    C200 = 1.2; C020 = 1.5; C002 = 1.8;
    C110 = 0.3; C101 = -0.2; C011 = 0.4;
    
    M_orig = InitializeM4_35(rho, u, v, w, C200, C110, C101, C020, C011, C002);
    
    m = moment_struct('from_vector', M_orig);
    C_array = M4toC4_3D(m.M000, m.M100, m.M200, m.M300, m.M400, ...
                        m.M010, m.M110, m.M210, m.M310, m.M020, m.M120, m.M220, ...
                        m.M030, m.M130, m.M040, m.M001, m.M101, m.M201, m.M301, ...
                        m.M002, m.M102, m.M202, m.M003, m.M103, m.M004, ...
                        m.M011, m.M111, m.M211, m.M021, m.M121, m.M031, ...
                        m.M012, m.M112, m.M013, m.M022);
    
    verifyEqual(testCase, size(C_array), [5 5 5], 'C_array should be 5x5x5');
    verifyTrue(testCase, all(isfinite(C_array(:))), 'All elements should be finite');
    
    fprintf('PASS: M4toC4_3D produces correct output test\n');
end

function test_standardized_moments_bounds(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    [~, S4] = M2CS4_35(M);
    
    S110 = S4(7);
    S101 = S4(17);
    S011 = S4(26);
    
    verifyLessThanOrEqual(testCase, abs(S110), 1.0 + testCase.TestData.tol, ...
                          'S110 should be in [-1, 1]');
    verifyLessThanOrEqual(testCase, abs(S101), 1.0 + testCase.TestData.tol, ...
                          'S101 should be in [-1, 1]');
    verifyLessThanOrEqual(testCase, abs(S011), 1.0 + testCase.TestData.tol, ...
                          'S011 should be in [-1, 1]');
    
    fprintf('PASS: Standardized moments bounds test\n');
end

function test_moment_idx(testCase)
    idx = moment_idx('M110');
    verifyEqual(testCase, idx, 7, 'M110 should be at index 7');
    
    idx = moment_idx('M200');
    verifyEqual(testCase, idx, 3, 'M200 should be at index 3');
    
    fprintf('PASS: moment_idx test\n');
end
