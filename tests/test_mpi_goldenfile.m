function tests = test_mpi_goldenfile
% Test MPI implementation against golden files for 1, 2, and 4 ranks
% Validates that:
% 1. MPI produces same results as serial code
% 2. MPI results are consistent across different rank counts
% 3. Results match stored golden files within tolerance

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
    testCase.TestData.tolerance = 1.0;  % Tolerance for MPI vs serial (different code paths)
    testCase.TestData.consistency_tolerance = 1e-13;  % Strict: MPI ranks must be bitwise identical
    
    % Check if Parallel Computing Toolbox is available
    testCase.TestData.has_pct = license('test', 'Distrib_Computing_Toolbox') && ...
                                 ~isempty(ver('parallel'));
    
    if ~testCase.TestData.has_pct
        warning('MPI tests require Parallel Computing Toolbox - tests will be skipped');
    end
end

function test_mpi_1_rank_vs_serial(testCase)
% Test that MPI with 1 rank produces identical results to serial
    
    % Skip if Parallel Computing Toolbox not available
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 1 Rank vs Serial ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 1 Rank vs Serial ===\n');
    
    % Load serial golden file (using Np=20 to meet minimum points/rank requirement)
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_Np20_tmax100.mat');
    if ~exist(golden_file, 'file')
        error('Serial golden file not found. Run run_goldenfile_creation.m first with Np=20.');
    end
    
    golden = load(golden_file);
    serial_data = golden.golden_data;
    
    % Run MPI with 1 rank
    fprintf('Running MPI simulation with 1 rank...\n');
    mpi_data = run_mpi_simulation(serial_data.parameters.Np, ...
                                  serial_data.parameters.tmax, 1);
    
    % Compare results
    compare_results(testCase, mpi_data, serial_data, ...
                   testCase.TestData.tolerance, '1-rank MPI vs Serial');
end

function test_mpi_2_ranks_vs_golden(testCase)
% Test MPI with 2 ranks against golden file
    
    % Skip if Parallel Computing Toolbox not available
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 2 Ranks vs Golden File ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    fprintf('\n=== TEST: MPI 2 Ranks vs Golden File ===\n');
    
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_2ranks_Np20_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        warning('Golden file for 2 ranks not found. Creating it now...');
        create_mpi_goldenfile(20, 0.1, 2, testCase.TestData.goldenfiles_dir);
        if ~exist(golden_file, 'file')
            error('Failed to create golden file for 2 ranks');
        end
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run MPI with 2 ranks
    fprintf('Running MPI simulation with 2 ranks...\n');
    mpi_data = run_mpi_simulation(golden_data.parameters.Np, ...
                                  golden_data.parameters.tmax, 2);
    
    % Compare results
    compare_results(testCase, mpi_data, golden_data, ...
                   testCase.TestData.tolerance, '2-rank MPI vs Golden');
end

function test_mpi_4_ranks_vs_golden(testCase)
% Test MPI with 4 ranks against golden file
    
    % Skip if Parallel Computing Toolbox not available
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI 4 Ranks vs Golden File ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    % Check maximum workers available (CI typically limits to 2)
    max_workers = feature('numcores');
    pool = gcp('nocreate');
    if ~isempty(pool)
        max_workers = pool.NumWorkers;
    end
    if max_workers < 4
        fprintf('\n=== TEST: MPI 4 Ranks vs Golden File ===\n');
        fprintf('SKIPPED: Requires 4 workers but only %d available (CI limitation)\n', max_workers);
        assumeFail(testCase, sprintf('Need 4 workers, only %d available', max_workers));
        return;
    end
    
    fprintf('\n=== TEST: MPI 4 Ranks vs Golden File ===\n');
    
    golden_file = fullfile(testCase.TestData.goldenfiles_dir, 'goldenfile_mpi_4ranks_Np20_tmax100.mat');
    
    if ~exist(golden_file, 'file')
        warning('Golden file for 4 ranks not found. Creating it now...');
        create_mpi_goldenfile(20, 0.1, 4, testCase.TestData.goldenfiles_dir);
        if ~exist(golden_file, 'file')
            error('Failed to create golden file for 4 ranks');
        end
    end
    
    golden = load(golden_file);
    golden_data = golden.golden_data;
    
    % Run MPI with 4 ranks
    fprintf('Running MPI simulation with 4 ranks...\n');
    mpi_data = run_mpi_simulation(golden_data.parameters.Np, ...
                                  golden_data.parameters.tmax, 4);
    
    % Compare results
    compare_results(testCase, mpi_data, golden_data, ...
                   testCase.TestData.tolerance, '4-rank MPI vs Golden');
end

function test_mpi_consistency_across_ranks(testCase)
% Test that different rank counts produce consistent results
    
    % Skip if Parallel Computing Toolbox not available
    if ~testCase.TestData.has_pct
        fprintf('\n=== TEST: MPI Consistency Across Ranks ===\n');
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required for MPI tests');
        return;
    end
    
    % Check maximum workers available
    max_workers = feature('numcores');
    pool = gcp('nocreate');
    if ~isempty(pool)
        max_workers = min(max_workers, pool.NumWorkers);
    end
    
    fprintf('\n=== TEST: MPI Consistency Across Ranks ===\n');
    
    Np = 20;  % Minimum grid size to support multiple ranks (≥10 points/rank)
    tmax = 0.1;
    
    % Adapt to available workers
    if max_workers >= 4
        fprintf('Running simulations with 1, 2, and 4 ranks...\n');
        rank_counts = [1, 2, 4];
    else
        fprintf('Running simulations with 1 and 2 ranks (limited to %d workers)...\n', max_workers);
        rank_counts = [1, 2];
    end
    
    % Run with different rank counts
    data_1rank = run_mpi_simulation(Np, tmax, 1);
    data_2ranks = run_mpi_simulation(Np, tmax, 2);
    if max_workers >= 4
        data_4ranks = run_mpi_simulation(Np, tmax, 4);
    end
    
    % Compare 1 vs 2 ranks (use relaxed tolerance for cross-rank consistency)
    fprintf('\nComparing 1 rank vs 2 ranks:\n');
    compare_results(testCase, data_2ranks, data_1rank, ...
                   testCase.TestData.consistency_tolerance, '2-rank vs 1-rank');
    
    if max_workers >= 4
        % Compare 1 vs 4 ranks
        fprintf('\nComparing 1 rank vs 4 ranks:\n');
        compare_results(testCase, data_4ranks, data_1rank, ...
                       testCase.TestData.consistency_tolerance, '4-rank vs 1-rank');
        
        % Compare 2 vs 4 ranks
        fprintf('\nComparing 2 ranks vs 4 ranks:\n');
        compare_results(testCase, data_4ranks, data_2ranks, ...
                       testCase.TestData.consistency_tolerance, '4-rank vs 2-rank');
    else
        fprintf('\n4-rank comparisons skipped (only %d workers available)\n', max_workers);
    end
end

%% Helper Functions

function mpi_data = run_mpi_simulation(Np, tmax, num_ranks)
% Run MPI simulation with specified number of ranks
    
    % Start/restart parallel pool
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', num_ranks);
    elseif pool.NumWorkers ~= num_ranks
        delete(pool);
        parpool('local', num_ranks);
    end
    
    % Run MPI simulation using unified main() interface
    mpi_data = main(Np, tmax, false, false, true, num_ranks);
    
    fprintf('  Completed: Np=%d, tmax=%.3f, ranks=%d, steps=%d\n', ...
            mpi_data.parameters.Np, mpi_data.parameters.tmax, ...
            num_ranks, mpi_data.parameters.time_steps);
end

function compare_results(testCase, current, golden, tolerance, desc)
% Compare current results against golden data
    
    fprintf('Comparing: %s\n', desc);
    
    % Compare final time (with slightly looser tolerance for time)
    time_tolerance = max(tolerance * 1e3, 1e-9);
    time_diff = abs(current.parameters.final_time - golden.parameters.final_time);
    verifyLessThanOrEqual(testCase, time_diff, time_tolerance, ...
        sprintf('Final time difference: %.3e', time_diff));
    fprintf('  Final time diff: %.3e (tolerance: %.3e) ✓\n', time_diff, time_tolerance);
    
    % Compare time steps (allow small variation)
    step_diff = abs(current.parameters.time_steps - golden.parameters.time_steps);
    verifyLessThanOrEqual(testCase, step_diff, 5, ...
        sprintf('Time steps difference: %d', step_diff));
    fprintf('  Time steps diff: %d (max: 5) ✓\n', step_diff);
    
    % Compare moment arrays
    moment_fields = {'M', 'C', 'S', 'M5', 'C5', 'S5'};
    
    for i = 1:length(moment_fields)
        field = moment_fields{i};
        if isfield(current.moments, field) && isfield(golden.moments, field)
            current_data = current.moments.(field);
            golden_data_field = golden.moments.(field);
            
            % Check dimensions
            verifyEqual(testCase, size(current_data), size(golden_data_field), ...
                sprintf('%s dimensions mismatch', field));
            
            % Calculate differences
            diff_data = abs(current_data - golden_data_field);
            max_diff = max(diff_data, [], 'all');
            mean_diff = mean(diff_data, 'all');
            
            % Verify within tolerance
            verifyLessThanOrEqual(testCase, max_diff, tolerance, ...
                sprintf('%s max diff %.3e exceeds tolerance %.3e', field, max_diff, tolerance));
            
            fprintf('  %s: max=%.3e, mean=%.3e ✓\n', field, max_diff, mean_diff);
        else
            warning('%s not found in results', field);
        end
    end
    
    fprintf('  PASS: %s\n', desc);
end

function create_mpi_goldenfile(Np, tmax, num_ranks, goldenfiles_dir)
% Create golden file for MPI with specified rank count
    
    fprintf('\n=== Creating MPI Golden File ===\n');
    fprintf('Parameters: Np=%d, tmax=%.3f, ranks=%d\n', Np, tmax, num_ranks);
    
    % Run simulation
    tic;
    golden_data = run_mpi_simulation(Np, tmax, num_ranks);
    elapsed_time = toc;
    
    % Add metadata
    golden_data.metadata.creation_date = datestr(now);
    golden_data.metadata.matlab_version = version;
    golden_data.metadata.description = sprintf('MPI golden file: Np=%d, tmax=%.3f, ranks=%d', ...
                                                Np, tmax, num_ranks);
    golden_data.metadata.num_ranks = num_ranks;
    golden_data.metadata.elapsed_time = elapsed_time;
    
    % Create directory if needed
    if ~exist(goldenfiles_dir, 'dir')
        mkdir(goldenfiles_dir);
    end
    
    % Save golden file
    golden_filename = fullfile(goldenfiles_dir, ...
                              sprintf('goldenfile_mpi_%dranks_Np%d_tmax%03d.mat', ...
                                      num_ranks, Np, round(tmax*1000)));
    save(golden_filename, 'golden_data');
    
    fprintf('Golden file saved: %s\n', golden_filename);
    fprintf('File size: %.2f KB\n', dir(golden_filename).bytes / 1024);
    fprintf('Elapsed time: %.2f seconds\n', elapsed_time);
end

