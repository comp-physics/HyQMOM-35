function tests = test_flux_eigenvalues
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-10;
end

function test_Flux_closure35_basic(testCase)
    rho = 1.0;
    u = 0.1; v = 0.2; w = 0.3;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [Fx, Fy, Fz, M_real] = Flux_closure35_and_realizable_3D(M, 1, 0.5);
    
    verifyEqual(testCase, length(Fx), 35, 'Fx should have 35 elements');
    verifyEqual(testCase, length(Fy), 35, 'Fy should have 35 elements');
    verifyEqual(testCase, length(Fz), 35, 'Fz should have 35 elements');
    
    verifyTrue(testCase, all(isfinite(Fx)), 'All Fx elements should be finite');
    verifyTrue(testCase, all(isfinite(Fy)), 'All Fy elements should be finite');
    verifyTrue(testCase, all(isfinite(Fz)), 'All Fz elements should be finite');
    
    fprintf('PASS: Flux_closure35 basic test\n');
end

function test_Flux_closure35_realizability(testCase)
    rho = 1.5;
    u = 0.0; v = 0.0; w = 0.0;
    T = 1.2;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [~, ~, ~, M_real] = Flux_closure35_and_realizable_3D(M, 1, 0.5);
    
    [realizable, violations] = verify_realizability(M_real, testCase.TestData.tol);
    
    verifyTrue(testCase, realizable, ...
               sprintf('Realized moments should be realizable. Violations: %d', violations.count));
    
    fprintf('PASS: Flux_closure35 realizability test\n');
end

function test_eigenvalues6_x_direction(testCase)
    rho = 1.0;
    u = 0.5; v = 0.0; w = 0.0;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [v6min, v6max, M_corr] = eigenvalues6_hyperbolic_3D(M, 1, 1, 0.5);
    
    verifyTrue(testCase, isreal(v6min), 'v6min should be real');
    verifyTrue(testCase, isreal(v6max), 'v6max should be real');
    verifyLessThanOrEqual(testCase, v6min, v6max, 'v6min <= v6max');
    
    fprintf('PASS: eigenvalues6 x-direction test\n');
end

function test_eigenvalues6_y_direction(testCase)
    rho = 1.0;
    u = 0.0; v = 0.5; w = 0.0;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [v6min, v6max, M_corr] = eigenvalues6_hyperbolic_3D(M, 2, 1, 0.5);
    
    verifyTrue(testCase, isreal(v6min), 'v6min should be real');
    verifyTrue(testCase, isreal(v6max), 'v6max should be real');
    verifyLessThanOrEqual(testCase, v6min, v6max, 'v6min <= v6max');
    
    fprintf('PASS: eigenvalues6 y-direction test\n');
end

function test_axis_moment_slice(testCase)
    rho = 1.0;
    u = 0.1; v = 0.2; w = 0.3;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    moments_uv = axis_moment_slice(M, 1);
    moments_vu = axis_moment_slice(M, 2);
    moments_uw = axis_moment_slice(M, 3);
    moments_vw = axis_moment_slice(M, 4);
    
    verifyEqual(testCase, length(moments_uv), 15, 'UV plane should have 10 moments');
    verifyEqual(testCase, length(moments_vu), 15, 'VU plane should have 10 moments');
    verifyEqual(testCase, length(moments_uw), 15, 'UW plane should have 10 moments');
    verifyEqual(testCase, length(moments_vw), 15, 'VW plane should have 10 moments');
    
    fprintf('PASS: axis_moment_slice test\n');
end

function test_flux_mass_flux_consistency(testCase)
    rho = 2.0;
    u = 0.3; v = 0.0; w = 0.0;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [Fx, ~, ~, ~] = Flux_closure35_and_realizable_3D(M, 1, 0.5);
    
    mass_flux_expected = rho * u;
    mass_flux_computed = Fx(1);
    
    verifyEqual(testCase, mass_flux_computed, mass_flux_expected, 'AbsTol', testCase.TestData.tol, ...
                'Mass flux should equal rho*u');
    
    fprintf('PASS: flux mass flux consistency test\n');
end
