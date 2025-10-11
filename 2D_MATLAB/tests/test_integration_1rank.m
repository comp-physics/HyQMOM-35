function tests = test_integration_1rank
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(script_dir);  % Add parent directory for main.m
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-6;
    testCase.TestData.has_pct = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
end

function test_small_grid_10x10(testCase)
    if ~testCase.TestData.has_pct
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required');
        return;
    end
    
    fprintf('\n=== Integration Test: 10x10 grid, 1 rank, 5 steps ===\n');
    
    Np = 10;
    tmax = 0.001;
    num_workers = 1;
    
    results = main(Np, tmax, false, num_workers, false, 1, false);
    
    verifyTrue(testCase, results.parameters.time_steps >= 1, 'Should take at least 1 time step');
    verifyEqual(testCase, size(results.moments.M, 1), Np, 'Grid size should match');
    verifyEqual(testCase, size(results.moments.M, 2), Np, 'Grid size should match');
    
    verifyTrue(testCase, all(results.moments.M(:,:,1) > 0, 'all'), 'Density should remain positive');
    
    [realizable, violations] = verify_realizability(results.moments.M, testCase.TestData.tol);
    verifyTrue(testCase, realizable, sprintf('Moments should be realizable. Violations: %d', violations.count));
    
    fprintf('PASS: Small grid 10x10 integration test\n');
end

function test_medium_grid_20x20(testCase)
    if ~testCase.TestData.has_pct
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required');
        return;
    end
    
    fprintf('\n=== Integration Test: 20x20 grid, 1 rank, 3 steps ===\n');
    
    Np = 20;
    tmax = 0.0005;
    num_workers = 1;
    
    results = main(Np, tmax, false, num_workers, false, 1, false);
    
    verifyTrue(testCase, results.parameters.time_steps >= 1, 'Should take at least 1 time step');
    verifyEqual(testCase, size(results.moments.M, 1), Np, 'Grid size should match');
    
    verifyTrue(testCase, all(results.moments.M(:,:,1) > 0, 'all'), 'Density should remain positive');
    
    fprintf('PASS: Medium grid 20x20 integration test\n');
end

function test_symmetry_preservation(testCase)
    if ~testCase.TestData.has_pct
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required');
        return;
    end
    
    fprintf('\n=== Integration Test: Symmetry preservation ===\n');
    
    Np = 10;
    tmax = 0.0005;
    num_workers = 1;
    
    results = main(Np, tmax, false, num_workers, false, 1, false);
    
    M = results.moments.M;
    
    diag_vals = zeros(Np, 5);
    for i = 1:Np
        diag_vals(i, 1) = M(i, i, 1);
        diag_vals(i, 2) = M(i, i, 2);
        diag_vals(i, 3) = M(i, i, 3);
        diag_vals(i, 4) = M(i, i, 4);
        diag_vals(i, 5) = M(i, i, 5);
    end
    
    Diff = zeros(Np, 5);
    for i = 1:Np
        j = Np + 1 - i;
        Diff(i, 1) = diag_vals(i, 1) - diag_vals(j, 1);
        Diff(i, 2) = diag_vals(i, 2) + diag_vals(j, 2);
        Diff(i, 3) = diag_vals(i, 3) - diag_vals(j, 3);
        Diff(i, 4) = diag_vals(i, 4) + diag_vals(j, 4);
        Diff(i, 5) = diag_vals(i, 5) - diag_vals(j, 5);
    end
    
    MaxDiff = zeros(5, 1);
    for k = 1:5
        Normk = norm(Diff(:, k));
        MaxDiff(k) = max(abs(Diff(:, k))) / (Normk + 1);
    end
    
    max_symmetry_error = max(abs(MaxDiff));
    fprintf('  Max symmetry error: %.6e\n', max_symmetry_error);
    
    verifyLessThan(testCase, max_symmetry_error, 0.1, 'Symmetry should be reasonably preserved');
    
    fprintf('PASS: Symmetry preservation test\n');
end
