function tests = test_mpi_local_only
%TEST_MPI_LOCAL_ONLY Local-only MPI tests for 4 and 8 ranks
%   These tests are SKIPPED in CI environments (detected by checking for
%   CI environment variables). They test larger MPI configurations that
%   require more workers than typically available in CI.
%   To run locally:
%     cd tests
%     runtests('test_mpi_local_only')
%   Note: Run create_goldenfiles_local.m first to generate golden files.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Setup for all tests
    
    % Add parent directory to path
    addpath('..');
    % Add src directory to path
    addpath('../src');
    addpath('../src/autogen');
    
    % Store paths in test case data
    testCase.TestData.goldenfiles_dir = '../goldenfiles';
    testCase.TestData.tolerance = 1e-6;  % Tolerance for MPI (floating-point differences)
    
    % Check if Parallel Computing Toolbox is available
    testCase.TestData.has_pct = license('test', 'Distrib_Computing_Toolbox') && ...
                                 ~isempty(ver('parallel'));
    
    % Detect CI environment
    testCase.TestData.is_ci = is_ci_environment();
    
    if testCase.TestData.is_ci
        fprintf('\n+==============================================================+\n');
        fprintf('|  CI ENVIRONMENT DETECTED - SKIPPING LOCAL-ONLY TESTS        |\n');
        fprintf('+==============================================================+\n\n');
    end
    
    if ~testCase.TestData.has_pct
        warning('MPI tests require Parallel Computing Toolbox - tests will be skipped');
    end
end

function test_mpi_4_ranks_vs_golden(testCase)
    % Test MPI with 4 ranks against golden file (40x40 grid)
    
    % Skip in CI
    if testCase.TestData.is_ci
        fprintf('SKIPPED: Local-only test (4 ranks not available in CI)\n');
        assumeFail(testCase, 'Local-only test - skipped in CI');
        return;
    end
    
    if ~testCase.TestData.has_pct
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required');
        return;
    end
    
    fprintf('\n=== TEST: MPI 4 Ranks (40x40 grid) [LOCAL ONLY] ===\n');
    
    num_ranks = 4;
    Np = 40;
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_4ranks_Np40_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        error(['Golden file for 4 ranks not found.\n' ...
               'Run create_goldenfiles_local.m first to generate it.']);
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run simulation
    fprintf('  Running simulation with %d ranks (Np=%d)...\n', num_ranks, Np);
    results = main(Np, golden_data.tmax, false, num_ranks, false);
    
    % Compare results
    fprintf('  Comparing against golden file...\n');
    
    % Extract final state
    if iscell(results.moments.M)
        M_final = results.moments.M{1};
    else
        M_final = results.moments.M;
    end
    
    % Compute max differences
    diff_M = max(abs(M_final(:) - golden_data.M(:)));
    
    fprintf('    max diff M: %.2e\n', diff_M);
    fprintf('    final time: %.6f (golden: %.6f)\n', results.parameters.final_time, golden_data.final_time);
    
    % Verify
    verifyLessThanOrEqual(testCase, diff_M, testCase.TestData.tolerance, ...
        sprintf('Moment differences exceed tolerance for %d ranks', num_ranks));
    
    fprintf('  OK 4-rank test PASSED\n');
end

function test_mpi_8_ranks_vs_golden(testCase)
    % Test MPI with 8 ranks against golden file (40x40 grid)
    
    % Skip in CI
    if testCase.TestData.is_ci
        fprintf('SKIPPED: Local-only test (8 ranks not available in CI)\n');
        assumeFail(testCase, 'Local-only test - skipped in CI');
        return;
    end
    
    if ~testCase.TestData.has_pct
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required');
        return;
    end
    
    fprintf('\n=== TEST: MPI 8 Ranks (40x40 grid) [LOCAL ONLY] ===\n');
    
    num_ranks = 8;
    Np = 40;
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_8ranks_Np40_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        error(['Golden file for 8 ranks not found.\n' ...
               'Run create_goldenfiles_local.m first to generate it.']);
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run simulation
    fprintf('  Running simulation with %d ranks (Np=%d)...\n', num_ranks, Np);
    results = main(Np, golden_data.tmax, false, num_ranks, false);
    
    % Compare results
    fprintf('  Comparing against golden file...\n');
    
    % Extract final state
    if iscell(results.moments.M)
        M_final = results.moments.M{1};
    else
        M_final = results.moments.M;
    end
    
    % Compute max differences
    diff_M = max(abs(M_final(:) - golden_data.M(:)));
    
    fprintf('    max diff M: %.2e\n', diff_M);
    fprintf('    final time: %.6f (golden: %.6f)\n', results.parameters.final_time, golden_data.final_time);
    
    % Verify
    verifyLessThanOrEqual(testCase, diff_M, testCase.TestData.tolerance, ...
        sprintf('Moment differences exceed tolerance for %d ranks', num_ranks));
    
    fprintf('  OK 8-rank test PASSED\n');
end

%% Helper function to detect CI environment
function is_ci = is_ci_environment()
    % Check common CI environment variables
    ci_vars = {'CI', 'CONTINUOUS_INTEGRATION', 'GITHUB_ACTIONS', ...
               'GITLAB_CI', 'CIRCLECI', 'TRAVIS', 'JENKINS_URL', ...
               'BUILDKITE', 'TEAMCITY_VERSION'};
    
    is_ci = false;
    for i = 1:length(ci_vars)
        if ~isempty(getenv(ci_vars{i}))
            is_ci = true;
            return;
        end
    end
end

