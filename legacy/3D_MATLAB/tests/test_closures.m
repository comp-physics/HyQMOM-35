function tests = test_closures
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-12;
end

function test_hyqmom_3D_gaussian(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.0;
    
    M4 = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    [~, S4] = M2CS4_35(M4);
    
    S300=S4(4);  S400=S4(5);  S110=S4(7);  S210=S4(8);  S310=S4(9);
    S120=S4(11); S220=S4(12); S030=S4(13); S130=S4(14); S040=S4(15);
    S101=S4(17); S201=S4(18); S301=S4(19); S102=S4(21); S202=S4(22);
    S003=S4(23); S103=S4(24); S004=S4(25); S011=S4(26); S111=S4(27);
    S211=S4(28); S021=S4(29); S121=S4(30); S031=S4(31); S012=S4(32);
    S112=S4(33); S013=S4(34); S022=S4(35);
    
    [S500,S410,S320,S230,S140,S401,S302,S203,S104,S311,S221,S131,S212,S113,S122,...
     S050,S041,S032,S023,S014,S005] = ...
        hyqmom_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                  S101,S201,S301,S102,S202,S003,S103,S004,...
                  S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
    
    verifyEqual(testCase, S500, 0.0, 'AbsTol', testCase.TestData.tol, 'S500 should be 0 for Gaussian');
    verifyEqual(testCase, S050, 0.0, 'AbsTol', testCase.TestData.tol, 'S050 should be 0 for Gaussian');
    verifyEqual(testCase, S005, 0.0, 'AbsTol', testCase.TestData.tol, 'S005 should be 0 for Gaussian');
    
    fprintf('PASS: hyqmom_3D Gaussian test\n');
end

function test_hyqmom_3D_all_closures_finite(testCase)
    rho = 1.5;
    u = 0.1; v = 0.2; w = 0.3;
    T = 1.2;
    
    M4 = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    [~, S4] = M2CS4_35(M4);
    
    S300=S4(4);  S400=S4(5);  S110=S4(7);  S210=S4(8);  S310=S4(9);
    S120=S4(11); S220=S4(12); S030=S4(13); S130=S4(14); S040=S4(15);
    S101=S4(17); S201=S4(18); S301=S4(19); S102=S4(21); S202=S4(22);
    S003=S4(23); S103=S4(24); S004=S4(25); S011=S4(26); S111=S4(27);
    S211=S4(28); S021=S4(29); S121=S4(30); S031=S4(31); S012=S4(32);
    S112=S4(33); S013=S4(34); S022=S4(35);
    
    [S500,S410,S320,S230,S140,S401,S302,S203,S104,S311,S221,S131,S212,S113,S122,...
     S050,S041,S032,S023,S014,S005] = ...
        hyqmom_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                  S101,S201,S301,S102,S202,S003,S103,S004,...
                  S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
    
    all_closures = [S500,S410,S320,S230,S140,S401,S302,S203,S104,S311,S221,S131,...
                    S212,S113,S122,S050,S041,S032,S023,S014,S005];
    
    verifyTrue(testCase, all(isfinite(all_closures)), 'All 21 closures should be finite');
    
    fprintf('PASS: hyqmom_3D all closures finite test\n');
end

function test_closure_and_eigenvalues_gaussian(testCase)
    mom = [1.0, 0.0, 1.0, 0.0, 3.0];
    
    [Mp, vpmin, vpmax] = closure_and_eigenvalues(mom);
    
    verifyTrue(testCase, isfinite(Mp), 'Closure moment should be finite');
    verifyTrue(testCase, isfinite(vpmin), 'Min eigenvalue should be finite');
    verifyTrue(testCase, isfinite(vpmax), 'Max eigenvalue should be finite');
    verifyLessThanOrEqual(testCase, vpmin, vpmax, 'vpmin should be less than vpmax');
    
    fprintf('PASS: closure_and_eigenvalues Gaussian test\n');
end

function test_Moments5_3D_consistency(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.0;
    
    M4 = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [M5, C5, S5] = Moments5_3D(M4);
    
    verifyEqual(testCase, length(M5), 21, 'M5 should have 21 elements');
    verifyEqual(testCase, length(C5), 21, 'C5 should have 21 elements');
    verifyEqual(testCase, length(S5), 21, 'S5 should have 21 elements');
    
    verifyTrue(testCase, all(isfinite(M5)), 'All M5 elements should be finite');
    verifyTrue(testCase, all(isfinite(C5)), 'All C5 elements should be finite');
    verifyTrue(testCase, all(isfinite(S5)), 'All S5 elements should be finite');
    
    fprintf('PASS: Moments5_3D consistency test\n');
end

function test_Moments5_3D_gaussian_symmetry(testCase)
    rho = 1.0;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.0;
    
    M4 = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    [~, ~, S5] = Moments5_3D(M4);
    
    verifyEqual(testCase, S5(1), 0.0, 'AbsTol', testCase.TestData.tol, 'S500 should be 0');
    verifyEqual(testCase, S5(16), 0.0, 'AbsTol', testCase.TestData.tol, 'S050 should be 0');
    verifyEqual(testCase, S5(21), 0.0, 'AbsTol', testCase.TestData.tol, 'S005 should be 0');
    
    fprintf('PASS: Moments5_3D Gaussian symmetry test\n');
end

function test_closure_eigenvalue_bounds(testCase)
    mom = [1.0, 0.5, 1.5, 0.3, 4.0];
    
    [~, vpmin, vpmax] = closure_and_eigenvalues(mom);
    
    verifyTrue(testCase, isreal(vpmin), 'vpmin should be real');
    verifyTrue(testCase, isreal(vpmax), 'vpmax should be real');
    verifyLessThanOrEqual(testCase, vpmin, vpmax, 'Eigenvalue ordering');
    
    fprintf('PASS: closure eigenvalue bounds test\n');
end
