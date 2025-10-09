function exit_code = run_all_tests()
%RUN_ALL_TESTS Master test runner for CI/CD
%   Runs all test files in sequence and reports results
%   Returns 0 if all tests pass, 1 if any fail

fprintf('\n');
fprintf('========================================\n');
fprintf('  HyQMOM-3D Comprehensive Test Suite\n');
fprintf('========================================\n');
fprintf('\n');

% Add paths
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'src'));
addpath(fullfile(script_dir, '..', 'src', 'autogen'));
addpath(fullfile(script_dir, 'utils'));

% List of all test files (in order)
test_files = {
    'test_moment_conversions';
    'test_initialization';
    'test_realizability';
    'test_closures';
    'test_flux_eigenvalues';
    'test_numerical_schemes';
    'test_collision';
    'test_diagnostics';
    'test_grid_utils';
    'test_autogen';
    'test_edge_corner_correction';
    'test_integration_1rank';
};

% Track results
total_tests = 0;
total_passed = 0;
total_failed = 0;
failed_tests = {};

start_time = tic;

% Run each test file
for i = 1:length(test_files)
    test_name = test_files{i};
    
    fprintf('\n----------------------------------------\n');
    fprintf('Running: %s\n', test_name);
    fprintf('----------------------------------------\n');
    
    try
        % Check if test file exists
        test_path = fullfile(script_dir, [test_name '.m']);
        if ~exist(test_path, 'file')
            fprintf('WARNING: Test file not found: %s\n', test_path);
            continue;
        end
        
        % Run the test
        result = runtests(test_name);
        
        % Count results
        num_passed = sum([result.Passed]);
        num_failed = sum([result.Failed]);
        num_incomplete = sum([result.Incomplete]);
        
        total_tests = total_tests + length(result);
        total_passed = total_passed + num_passed;
        total_failed = total_failed + num_failed + num_incomplete;
        
        % Report results for this file
        if num_failed == 0 && num_incomplete == 0
            fprintf('OK PASSED: %d/%d tests\n', num_passed, length(result));
        else
            fprintf('X FAILED: %d passed, %d failed, %d incomplete\n', ...
                    num_passed, num_failed, num_incomplete);
            failed_tests{end+1} = test_name;
            
            % Show details of failures
            for j = 1:length(result)
                if result(j).Failed || result(j).Incomplete
                    fprintf('  - %s: %s\n', result(j).Name, result(j).Details.DiagnosticRecord.Report);
                end
            end
        end
        
    catch ME
        fprintf('X ERROR running %s:\n', test_name);
        fprintf('  %s\n', ME.message);
        failed_tests{end+1} = test_name;
        total_failed = total_failed + 1;
    end
end

elapsed_time = toc(start_time);

% Print summary
fprintf('\n');
fprintf('========================================\n');
fprintf('  Test Suite Summary\n');
fprintf('========================================\n');
fprintf('Total tests run:  %d\n', total_tests);
fprintf('Passed:           %d\n', total_passed);
fprintf('Failed:           %d\n', total_failed);
fprintf('Time elapsed:     %.2f seconds\n', elapsed_time);
fprintf('\n');

if total_failed == 0
    fprintf('OKOKOK ALL TESTS PASSED OKOKOK\n');
    exit_code = 0;
else
    fprintf('XXX SOME TESTS FAILED XXX\n');
    fprintf('\nFailed test files:\n');
    for i = 1:length(failed_tests)
        fprintf('  - %s\n', failed_tests{i});
    end
    exit_code = 1;
end

fprintf('========================================\n');
fprintf('\n');

% Exit with appropriate code for CI
if nargout == 0
    exit(exit_code);
end

end
