function tests = test_numerical_schemes
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-10;
end

function test_flux_HLL_basic(testCase)
    N = 10;
    Nmom = 5;
    
    Wstar = rand(N, Nmom);
    W = rand(N, Nmom);
    l1 = -ones(N, 1);
    l2 = ones(N, 1);
    F = rand(N, Nmom);
    
    Flux = flux_HLL(Wstar, W, l1, l2, F, N);
    
    verifyEqual(testCase, size(Flux, 1), N-1, 'Flux should have N-1 rows');
    verifyEqual(testCase, size(Flux, 2), Nmom, 'Flux should have Nmom columns');
    verifyTrue(testCase, all(isfinite(Flux(:))), 'All flux values should be finite');
    
    fprintf('PASS: flux_HLL basic test\n');
end

function test_pas_HLL_interior(testCase)
    Np = 20;
    Nmom = 35;
    
    M = ones(Np, Nmom);
    F = ones(Np, Nmom);
    dt = 0.01;
    dx = 0.1;
    vpmin = -ones(Np, 1);
    vpmax = ones(Np, 1);
    
    Mp = pas_HLL(M, F, dt, dx, vpmin, vpmax, true, true);
    
    verifyEqual(testCase, size(Mp), size(M), 'Output size should match input');
    verifyTrue(testCase, all(isfinite(Mp(:))), 'All updated moments should be finite');
    
    fprintf('PASS: pas_HLL interior test\n');
end

function test_pas_HLL_no_left_bc(testCase)
    Np = 20;
    Nmom = 35;
    
    M = ones(Np, Nmom);
    F = ones(Np, Nmom);
    dt = 0.01;
    dx = 0.1;
    vpmin = -ones(Np, 1);
    vpmax = ones(Np, 1);
    
    Mp = pas_HLL(M, F, dt, dx, vpmin, vpmax, false, true);
    
    verifyEqual(testCase, size(Mp), size(M), 'Output size should match input');
    verifyTrue(testCase, all(isfinite(Mp(:))), 'All updated moments should be finite');
    
    fprintf('PASS: pas_HLL no left BC test\n');
end

function test_pas_HLL_no_right_bc(testCase)
    Np = 20;
    Nmom = 35;
    
    M = ones(Np, Nmom);
    F = ones(Np, Nmom);
    dt = 0.01;
    dx = 0.1;
    vpmin = -ones(Np, 1);
    vpmax = ones(Np, 1);
    
    Mp = pas_HLL(M, F, dt, dx, vpmin, vpmax, true, false);
    
    verifyEqual(testCase, size(Mp), size(M), 'Output size should match input');
    verifyTrue(testCase, all(isfinite(Mp(:))), 'All updated moments should be finite');
    
    fprintf('PASS: pas_HLL no right BC test\n');
end

function test_pas_HLL_processor_boundaries(testCase)
    Np = 20;
    Nmom = 35;
    
    M = ones(Np, Nmom);
    F = ones(Np, Nmom);
    dt = 0.01;
    dx = 0.1;
    vpmin = -ones(Np, 1);
    vpmax = ones(Np, 1);
    
    Mp = pas_HLL(M, F, dt, dx, vpmin, vpmax, false, false);
    
    verifyEqual(testCase, size(Mp), size(M), 'Output size should match input');
    verifyTrue(testCase, all(isfinite(Mp(:))), 'All updated moments should be finite');
    
    fprintf('PASS: pas_HLL processor boundaries test\n');
end

function test_pas_HLL_stability(testCase)
    Np = 20;
    Nmom = 35;
    
    rho = 1.0;
    u = 0.1; v = 0.0; w = 0.0;
    T = 1.0;
    M_vec = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    M = repmat(M_vec', Np, 1);
    
    [Fx, ~, ~, ~] = Flux_closure35_and_realizable_3D(M_vec, 1, 0.5);
    F = repmat(Fx', Np, 1);
    
    dt = 0.001;
    dx = 0.1;
    vpmin = -2*ones(Np, 1);
    vpmax = 2*ones(Np, 1);
    
    Mp = pas_HLL(M, F, dt, dx, vpmin, vpmax, true, true);
    
    verifyTrue(testCase, all(Mp(:, 1) > 0), 'Density should remain positive');
    
    fprintf('PASS: pas_HLL stability test\n');
end
