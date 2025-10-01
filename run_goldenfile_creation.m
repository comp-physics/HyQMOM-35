% Script to create golden file for regression testing
% This creates the rigorous test with Np=10, tmax=0.1

fprintf('=== GOLDEN FILE CREATION ===\n');
fprintf('This script will create a golden file by running the simulation with:\n');
fprintf('  - Np = 10 (grid points)\n');
fprintf('  - tmax = 0.1 (final time)\n');
fprintf('  - enable_plots = false (no plotting)\n\n');
fprintf('This is a rigorous test that validates longer time integration.\n\n');

% Add src directory to path for function dependencies
addpath('src');

% Check if we're in the right directory
required_files = {'main.m', 'simulation_plots.m'};
for i = 1:length(required_files)
    if ~exist(required_files{i}, 'file')
        error('%s not found in current directory', required_files{i});
    end
end

% Create goldenfiles directory if it doesn't exist
goldenfiles_dir = 'goldenfiles';
if ~exist(goldenfiles_dir, 'dir')
    mkdir(goldenfiles_dir);
    fprintf('Created directory: %s\n', goldenfiles_dir);
end

% Define golden file parameters
GOLDEN_NP = 10;
GOLDEN_TMAX = 0.1;
GOLDEN_ENABLE_PLOTS = false;
GOLDEN_SAVE_OUTPUT = false;  % We manually save the golden data structure

% Run the golden file creation
try
    fprintf('Golden file parameters: Np = %d, tmax = %.3f, enable_plots = %s\n', ...
        GOLDEN_NP, GOLDEN_TMAX, mat2str(GOLDEN_ENABLE_PLOTS));
    fprintf('Running simulation...\n');
    tic;
    
    % Execute the main simulation
    results = main(GOLDEN_NP, GOLDEN_TMAX, GOLDEN_ENABLE_PLOTS, GOLDEN_SAVE_OUTPUT);
    
    elapsed_time = toc;
    fprintf('Simulation completed in %.2f seconds\n', elapsed_time);
    
    % Create golden file structure
    golden_data = results;
    golden_data.parameters.elapsed_time = elapsed_time;
    
    % Add metadata
    golden_data.metadata.creation_date = datestr(now);
    golden_data.metadata.matlab_version = version;
    golden_data.metadata.description = sprintf('Rigorous golden file: Np=%d, tmax=%.3f', GOLDEN_NP, GOLDEN_TMAX);
    golden_data.metadata.git_branch = 'refac-clean-baseline';
    
    % Save golden file
    golden_filename = fullfile(goldenfiles_dir, sprintf('goldenfile_Np%d_tmax%03d.mat', GOLDEN_NP, round(GOLDEN_TMAX*1000)));
    save(golden_filename, 'golden_data');
    
    fprintf('\n=== SUCCESS ===\n');
    fprintf('Golden file saved: %s\n', golden_filename);
    fprintf('File size: %.2f KB\n', dir(golden_filename).bytes / 1024);
    fprintf('Final time reached: %.4f\n', golden_data.parameters.final_time);
    fprintf('Time steps: %d\n', golden_data.parameters.time_steps);
    
    % Verify symmetry
    fprintf('\nVerifying symmetry in saved data:\n');
    M = golden_data.moments.M;
    for mom = 1:5
        maxdiff = max(abs(M(:,:,mom) - flipud(M(:,:,mom))), [], 'all');
        fprintf('  M%d max asymmetry: %.2e\n', mom-1, maxdiff);
    end
    
    fprintf('\nUse tests/test_validate_goldenfile.m to validate future changes.\n');
    
catch ME
    fprintf('\n=== ERROR ===\n');
    fprintf('ERROR during simulation: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
        end
    end
    fprintf('\nScript completed.\n');
    rethrow(ME);
end
