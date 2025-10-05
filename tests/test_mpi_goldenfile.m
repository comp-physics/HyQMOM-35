function tests = test_mpi_goldenfile
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
    
    if ~testCase.TestData.has_pct
        warning('MPI tests require Parallel Computing Toolbox - tests will be skipped');
    end
end

function test_mpi_1_rank_vs_golden(testCase)
    % Test MPI with 1 rank (single processor) against golden file (20×20 grid)
    
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 1 Rank vs Golden ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 1 Rank (20×20 grid) ===\n');
    
    num_ranks = 1;
    Np = 20;
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_1ranks_Np20_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        error('Golden file for 1 rank not found. Run create_goldenfiles.m first.');
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run MPI simulation
    fprintf('Running MPI simulation with 1 rank...\n');
    mpi_data = run_mpi_simulation(Np, golden_data.parameters.tmax, num_ranks);
    
    % Compare against golden file
    compare_results(testCase, mpi_data, golden_data, ...
                   testCase.TestData.tolerance, '1-rank MPI vs Golden');
end

function test_mpi_2_ranks_vs_golden(testCase)
% Test MPI with 2 ranks against golden file (40×40 grid)
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 2 Ranks vs Golden ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 2 Ranks (40×40 grid) ===\n');
    
    num_ranks = 2;
    Np = 40;
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_2ranks_Np40_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        error('Golden file for 2 ranks not found. Run create_goldenfiles.m first.');
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run MPI simulation
    fprintf('Running MPI simulation with 2 ranks...\n');
    mpi_data = run_mpi_simulation(Np, golden_data.parameters.tmax, num_ranks);
    
    % Compare against golden file
    compare_results(testCase, mpi_data, golden_data, ...
                   testCase.TestData.tolerance, '2-rank MPI vs Golden');
end

%% Helper Functions

function mpi_data = run_mpi_simulation(Np, tmax, num_ranks)
    mpi_data = main(Np, tmax, false, num_ranks);
end

function compare_results(testCase, data1, data2, tolerance, description)
    % Compare two result structures
    
    fprintf('\nComparing results: %s\n', description);
    
    % Extract moment arrays
    M1 = data1.moments.M;
    M2 = data2.moments.M;
    
    % Check sizes match
    if ~isequal(size(M1), size(M2))
        error('Result sizes do not match: %s vs %s', mat2str(size(M1)), mat2str(size(M2)));
    end
    
    % Compute differences
    diff_M = abs(M1 - M2);
    max_diff = max(diff_M(:));
    mean_diff = mean(diff_M(:));
    rel_diff = max_diff / max(abs(M2(:)));
    
    fprintf('  Max absolute difference: %.6e\n', max_diff);
    fprintf('  Mean absolute difference: %.6e\n', mean_diff);
    fprintf('  Max relative difference: %.6e\n', rel_diff);
    
    % Sample values at center
    [nx, ny, ~] = size(M1);
    ic = round(nx/2);
    jc = round(ny/2);
    fprintf('  Sample M(center,center,1): %.10e (data1) vs %.10e (data2)\n', ...
            M1(ic,jc,1), M2(ic,jc,1));
    
    % Assess result
    if max_diff < tolerance
        if max_diff < 1e-14
            fprintf('  PERFECT: Bitwise identical!\n');
        else
            fprintf('  PASS: Within tolerance (%.1e)\n', tolerance);
        end
    else
        fprintf('  FAIL: Exceeds tolerance (%.1e)\n', tolerance);
        error('Results differ by %.6e, tolerance %.6e', max_diff, tolerance);
    end
    
    % MATLAB unit test assertion
    verifyLessThan(testCase, max_diff, tolerance, ...
                   sprintf('%s: max difference %.6e exceeds tolerance %.6e', ...
                           description, max_diff, tolerance));
end
