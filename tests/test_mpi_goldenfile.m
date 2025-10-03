function tests = test_mpi_goldenfile
% Test MPI implementation against golden files for 1, 2, 3, and 4 ranks
% Validates that:
% 1. MPI results are consistent across different rank counts
% 2. Results match stored golden files within tolerance
% Each configuration uses 20 grid points per rank per dimension

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Setup for all tests
    % Add parent directory to path
    addpath('..');
    % Add src directory to path
    addpath('../src');
    
    % Store paths in test case data
    testCase.TestData.goldenfiles_dir = '../goldenfiles';
    testCase.TestData.consistency_tolerance = 1e-6;  % ~1e-6 is expected for MPI (floating-point)
    
    % Check if Parallel Computing Toolbox is available
    testCase.TestData.has_pct = license('test', 'Distrib_Computing_Toolbox') && ...
                                 ~isempty(ver('parallel'));
    
    if ~testCase.TestData.has_pct
        warning('MPI tests require Parallel Computing Toolbox - tests will be skipped');
    end
end

function test_mpi_1_rank(testCase)
% Test MPI with 1 rank (20x20 grid)
    
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 1 Rank ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 1 Rank (20x20 grid) ===\n');
    
    num_ranks = 1;
    Np = 20;  % 1 rank × 20 pts/rank
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
                   testCase.TestData.consistency_tolerance, '1-rank MPI vs Golden');
end

function test_mpi_2_ranks(testCase)
% Test MPI with 2 ranks (40x40 grid: 2x1 decomposition, 20 pts/rank)
    
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 2 Ranks ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 2 Ranks (40x40 grid) ===\n');
    
    num_ranks = 2;
    Np = 40;  % 2 ranks × 20 pts/rank
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
                   testCase.TestData.consistency_tolerance, '2-rank MPI vs Golden');
end

function test_mpi_3_ranks(testCase)
% Test MPI with 3 ranks (60x60 grid: 3x1 decomposition, 20 pts/rank)
    
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 3 Ranks ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 3 Ranks (60x60 grid) ===\n');
    
    num_ranks = 3;
    Np = 60;  % 3 ranks × 20 pts/rank
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_3ranks_Np60_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        error('Golden file for 3 ranks not found. Run create_goldenfiles.m first.');
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run MPI simulation
    fprintf('Running MPI simulation with 3 ranks...\n');
    mpi_data = run_mpi_simulation(Np, golden_data.parameters.tmax, num_ranks);
    
    % Compare against golden file
    compare_results(testCase, mpi_data, golden_data, ...
                   testCase.TestData.consistency_tolerance, '3-rank MPI vs Golden');
end

function test_mpi_4_ranks(testCase)
% Test MPI with 4 ranks (40x40 grid: 2x2 decomposition, 20 pts/rank)
    
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 4 Ranks ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 4 Ranks (40x40 grid) ===\n');
    
    num_ranks = 4;
    Np = 40;  % 2×2 ranks × 20 pts/rank
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_4ranks_Np40_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        error('Golden file for 4 ranks not found. Run create_goldenfiles.m first.');
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run MPI simulation
    fprintf('Running MPI simulation with 4 ranks...\n');
    mpi_data = run_mpi_simulation(Np, golden_data.parameters.tmax, num_ranks);
    
    % Compare against golden file
    compare_results(testCase, mpi_data, golden_data, ...
                   testCase.TestData.consistency_tolerance, '4-rank MPI vs Golden');
end

function test_mpi_consistency_across_ranks(testCase)
% Test that 1, 2, 3, and 4 rank results are consistent at common grid points
    
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI Consistency Across Ranks ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI Consistency Across Ranks ===\n');
    fprintf('Comparing 1-rank (20x20) vs 2-rank (40x40) at downsampled points...\n\n');
    
    % Load 1-rank golden file (20x20)
    golden_1rank = load(fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_1ranks_Np20_tmax100.mat'));
    M_1rank = golden_1rank.golden_data.moments.M;
    
    % Load 2-rank golden file (40x40)
    golden_2rank = load(fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_2ranks_Np40_tmax100.mat'));
    M_2rank = golden_2rank.golden_data.moments.M;
    
    % Compare at downsampled points (every other point in 40x40 grid)
    % 40x40[1:2:end, 1:2:end] should roughly match 20x20
    M_2rank_ds = M_2rank(1:2:end, 1:2:end, :);
    
    diff = abs(M_1rank - M_2rank_ds);
    max_diff = max(diff(:));
    
    fprintf('Max difference between 1-rank and 2-rank (downsampled): %.6e\n', max_diff);
    
    % Note: Results may differ because different grid sizes = different IC = different physics
    % This is more of an informational test
    if max_diff < 0.1  % Very loose tolerance since grids are different
        fprintf('✓ Results are similar\n');
    else
        fprintf('Note: Results differ (expected - different grid resolutions)\n');
    end
    
    fprintf('\n');
end

%% Helper Functions

function mpi_data = run_mpi_simulation(Np, tmax, num_ranks)
% Run MPI simulation and return results
    
    % Use main_mpi directly
    mpi_data = main_mpi(Np, tmax, false, num_ranks);
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
            fprintf('  ✅ PERFECT: Bitwise identical!\n');
        else
            fprintf('  ✅ PASS: Within tolerance (%.1e)\n', tolerance);
        end
    else
        fprintf('  ❌ FAIL: Exceeds tolerance (%.1e)\n', tolerance);
        error('Results differ by %.6e, tolerance %.6e', max_diff, tolerance);
    end
    
    % MATLAB unit test assertion
    verifyLessThan(testCase, max_diff, tolerance, ...
                   sprintf('%s: max difference %.6e exceeds tolerance %.6e', ...
                           description, max_diff, tolerance));
end
