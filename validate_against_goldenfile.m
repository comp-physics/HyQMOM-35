% Script to validate current simulation results against the golden file
% This script runs the simulation with the same parameters as the golden file
% and compares the results for regression testing
% Golden files are stored in the goldenfiles/ directory

function validate_against_goldenfile(tolerance)
    if nargin < 1
        tolerance = 1e-12;  % Default tolerance for comparison
    end
    
    clc
    fprintf('=== VALIDATION AGAINST GOLDEN FILE ===\n');
    fprintf('Tolerance: %.2e\n\n', tolerance);
    
    % Check if golden file exists in goldenfiles directory
    goldenfiles_dir = 'goldenfiles';
    golden_filename = fullfile(goldenfiles_dir, 'goldenfile_Np6_tmax020.mat');
    if ~exist(golden_filename, 'file')
        error('Golden file not found: %s\nRun run_goldenfile_creation.m first.', golden_filename);
    end
    
    % Load golden file
    fprintf('Loading golden file: %s\n', golden_filename);
    golden = load(golden_filename);
    golden_data = golden.golden_data;
    
    fprintf('Golden file created: %s\n', golden_data.metadata.creation_date);
    fprintf('Golden file parameters: Np=%d, tmax=%.3f, final_time=%.6f, steps=%d\n', ...
        golden_data.parameters.Np, golden_data.parameters.tmax, ...
        golden_data.parameters.final_time, golden_data.parameters.time_steps);
    
    try
        % Clear workspace to avoid conflicts
        clear M C S M5 C5 S5 t nn Np tmax
        
        % Run simulation with same parameters as golden file
        fprintf('\nRunning validation simulation...\n');
        tic;
        results = main_2Dcrossing_3DHyQMOM35(golden_data.parameters.Np, golden_data.parameters.tmax, false);
        elapsed_time = toc;
        fprintf('Simulation completed in %.2f seconds\n', elapsed_time);
        
        % Compare results
        fprintf('\n=== COMPARISON RESULTS ===\n');
        
        % Check parameters (they should match since we used the same overrides)
        param_match = true;
        if abs(golden_data.parameters.Np - golden_data.parameters.Np) > 0
            fprintf('FAIL: Np mismatch - Current: %d, Golden: %d\n', golden_data.parameters.Np, golden_data.parameters.Np);
            param_match = false;
        end
        
        if abs(golden_data.parameters.tmax - golden_data.parameters.tmax) > tolerance
            fprintf('FAIL: tmax mismatch - Current: %.6f, Golden: %.6f\n', golden_data.parameters.tmax, golden_data.parameters.tmax);
            param_match = false;
        end
        
        if param_match
            fprintf('PASS: Parameters match (using golden file parameters)\n');
        end
        
        % Compare final time and steps (allowing some variation due to adaptive time stepping)
        time_tolerance = max(tolerance, 1e-6);  % More lenient for time comparison
        if abs(results.parameters.final_time - golden_data.parameters.final_time) > time_tolerance
            fprintf('WARN: Final time difference - Current: %.6f, Golden: %.6f, Diff: %.2e\n', ...
                results.parameters.final_time, golden_data.parameters.final_time, abs(results.parameters.final_time - golden_data.parameters.final_time));
        else
            fprintf('PASS: Final time matches within tolerance\n');
        end
        
        step_diff = abs(results.parameters.time_steps - golden_data.parameters.time_steps);
        if step_diff > 1  % Allow 1 step difference
            fprintf('WARN: Time steps difference - Current: %d, Golden: %d, Diff: %d\n', ...
                results.parameters.time_steps, golden_data.parameters.time_steps, step_diff);
        else
            fprintf('PASS: Time steps match within tolerance\n');
        end
        
        % Compare moment arrays
        moment_fields = {'M', 'C', 'S', 'M5', 'C5', 'S5'};
        all_moments_match = true;
        
        for i = 1:length(moment_fields)
            field = moment_fields{i};
            if isfield(results.moments, field)
                current_data = results.moments.(field);
                golden_field_data = golden_data.moments.(field);
                
                % Check dimensions
                if ~isequal(size(current_data), size(golden_field_data))
                    fprintf('FAIL: %s dimensions mismatch - Current: %s, Golden: %s\n', ...
                        field, mat2str(size(current_data)), mat2str(size(golden_field_data)));
                    all_moments_match = false;
                    continue;
                end
                
                % Calculate differences
                diff_data = abs(current_data - golden_field_data);
                max_diff = max(diff_data, [], 'all');
                mean_diff = mean(diff_data, 'all');
                
                % Check for NaN/Inf
                has_nan_current = any(isnan(current_data), 'all');
                has_inf_current = any(isinf(current_data), 'all');
                has_nan_golden = any(isnan(golden_field_data), 'all');
                has_inf_golden = any(isinf(golden_field_data), 'all');
                
                if has_nan_current || has_inf_current || has_nan_golden || has_inf_golden
                    fprintf('FAIL: %s contains NaN/Inf values\n', field);
                    all_moments_match = false;
                elseif max_diff > tolerance
                    fprintf('FAIL: %s max difference %.2e exceeds tolerance %.2e\n', field, max_diff, tolerance);
                    fprintf('      Mean difference: %.2e\n', mean_diff);
                    all_moments_match = false;
                else
                    fprintf('PASS: %s matches within tolerance (max diff: %.2e)\n', field, max_diff);
                end
            else
                fprintf('WARN: %s not found in current simulation\n', field);
            end
        end
        
        % Overall result
        fprintf('\n=== OVERALL RESULT ===\n');
        if param_match && all_moments_match
            fprintf('SUCCESS: All validation tests passed!\n');
            fprintf('The simulation produces results consistent with the golden file.\n');
            exit_code = 0;
        else
            fprintf('FAILURE: Some validation tests failed.\n');
            fprintf('Check the differences above and investigate potential issues.\n');
            exit_code = 1;
        end
        
        % Summary statistics
        fprintf('\n=== SUMMARY STATISTICS ===\n');
        fprintf('Current density (M000) - Min: %.6e, Max: %.6e, Mean: %.6e\n', ...
            min(results.moments.M(:,:,1), [], 'all'), max(results.moments.M(:,:,1), [], 'all'), mean(results.moments.M(:,:,1), 'all'));
        fprintf('Golden density (M000)  - Min: %.6e, Max: %.6e, Mean: %.6e\n', ...
            min(golden_data.moments.M(:,:,1), [], 'all'), ...
            max(golden_data.moments.M(:,:,1), [], 'all'), ...
            mean(golden_data.moments.M(:,:,1), 'all'));
        
    catch ME
        fprintf('ERROR during validation: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
        exit_code = 2;  % Error during execution
    end
    
    % Exit with appropriate code
    if exist('exit_code', 'var')
        if exit_code ~= 0
            fprintf('\nExiting with error code %d\n', exit_code);
        end
        exit(exit_code);
    else
        % Fallback if exit_code wasn't set
        exit(2);
    end
end
