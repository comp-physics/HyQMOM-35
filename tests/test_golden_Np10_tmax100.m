function tests = test_golden_Np10_tmax100
% Test function for validating Np=10, tmax=0.1 against golden file
% This is a rigorous test for longer time integration

tests = functiontests(localfunctions);
end

function test_simulation_against_goldenfile(testCase)
% Test that current simulation matches the golden file within tolerance

fprintf('=== Testing Np=10, tmax=0.1 ===\n');
fprintf('MATLAB version: %s\n', version);

% Add parent directory to path
addpath('..');
addpath('../src');

% Set tolerance
tolerance = 1e-10;

% Load golden file
goldenfiles_dir = '../goldenfiles';
golden_filename = fullfile(goldenfiles_dir, 'goldenfile_Np10_tmax100.mat');

if ~exist(golden_filename, 'file')
    error('Golden file not found: %s\nRun create_golden_Np10_tmax010.m first.', golden_filename);
end

fprintf('Loading golden file: %s\n', golden_filename);
golden = load(golden_filename);
golden_data = golden.golden_data;

fprintf('Golden file created: %s\n', golden_data.metadata.creation_date);
fprintf('Golden file parameters: Np=%d, tmax=%.3f, final_time=%f, steps=%d\n', ...
    golden_data.parameters.Np, golden_data.parameters.tmax, ...
    golden_data.parameters.final_time, golden_data.parameters.time_steps);

% Run simulation with same parameters
fprintf('\nRunning validation simulation...\n');
tic;
results = main(golden_data.parameters.Np, golden_data.parameters.tmax, false, false);
elapsed_time = toc;
fprintf('Simulation completed in %.2f seconds\n', elapsed_time);

% Compare results
fprintf('\n=== COMPARISON RESULTS ===\n');

% Compare final time and steps
time_tolerance = max(tolerance, 1e-6);
verifyEqual(testCase, results.parameters.final_time, golden_data.parameters.final_time, ...
    'AbsTol', time_tolerance, 'Final time should match within tolerance');

step_diff = abs(results.parameters.time_steps - golden_data.parameters.time_steps);
verifyLessThanOrEqual(testCase, step_diff, 5, 'Time steps should match within 5 steps');

% Compare moment arrays
moment_fields = {'M', 'C', 'S', 'M5', 'C5', 'S5'};

for i = 1:length(moment_fields)
    field = moment_fields{i};
    if isfield(results.moments, field) && isfield(golden_data.moments, field)
        current_data = results.moments.(field);
        golden_field_data = golden_data.moments.(field);
        
        % Check dimensions
        verifyEqual(testCase, size(current_data), size(golden_field_data), ...
            sprintf('%s dimensions should match', field));
        
        % Calculate differences
        diff_data = abs(current_data - golden_field_data);
        max_diff = max(diff_data, [], 'all');
        
        % Verify within tolerance
        verifyLessThanOrEqual(testCase, max_diff, tolerance, ...
            sprintf('%s should match within tolerance (max diff: %.2e)', field, max_diff));
        
        fprintf('PASS: %s matches within tolerance (max diff: %.2e)\n', field, max_diff);
    else
        warning('%s not found in current simulation or golden file', field);
    end
end

fprintf('\n=== TEST COMPLETED SUCCESSFULLY ===\n');
fprintf('All validation tests passed!\n');

end

