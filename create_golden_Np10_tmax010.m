% Script to create golden file for Np=10, tmax=0.1
% This creates a rigorous test for longer time integration

fprintf('=== GOLDEN FILE CREATION: Np=10, tmax=0.1 ===\n');
fprintf('This will create a golden file for rigorous testing with:\n');
fprintf('  - Np = 10 (grid points)\n');
fprintf('  - tmax = 0.1 (final time)\n');
fprintf('  - enable_plots = false (no plotting)\n\n');

% Add src directory to path
addpath('src');
addpath('src/autogen');

% Create goldenfiles directory if needed
goldenfiles_dir = 'goldenfiles';
if ~exist(goldenfiles_dir, 'dir')
    mkdir(goldenfiles_dir);
end

% Define parameters
GOLDEN_NP = 10;
GOLDEN_TMAX = 0.1;

% Define golden filename
golden_filename = fullfile(goldenfiles_dir, 'goldenfile_Np10_tmax100.mat');

try
    fprintf('Running simulation with Np=%d, tmax=%.3f...\n', GOLDEN_NP, GOLDEN_TMAX);
    tic;
    
    % Run simulation
    results = main(GOLDEN_NP, GOLDEN_TMAX, false, false);
    
    elapsed_time = toc;
    fprintf('Simulation completed in %.2f seconds\n', elapsed_time);
    
    % Create golden data structure
    golden_data = results;
    golden_data.parameters.elapsed_time = elapsed_time;
    
    % Add metadata
    golden_data.metadata.creation_date = datestr(now);
    golden_data.metadata.matlab_version = version;
    golden_data.metadata.description = sprintf('Golden file for Np=%d, tmax=%.3f (rigorous test)', GOLDEN_NP, GOLDEN_TMAX);
    golden_data.metadata.git_commit = 'refac-clean-baseline';
    
    % Save
    save(golden_filename, 'golden_data');
    
    fprintf('\n=== SUCCESS ===\n');
    fprintf('Golden file saved: %s\n', golden_filename);
    fprintf('File size: %.2f KB\n', dir(golden_filename).bytes / 1024);
    fprintf('Final time reached: %.4f\n', golden_data.parameters.final_time);
    fprintf('Time steps: %d\n', golden_data.parameters.time_steps);
    
    % Verify symmetry
    fprintf('\nVerifying symmetry in saved data:\n');
    M = golden_data.moments.M;
    Np = size(M, 1);
    for mom = 1:5
        maxdiff = max(abs(M(:,:,mom) - flipud(M(:,:,mom))), [], 'all');
        fprintf('  M%d max asymmetry: %.2e\n', mom-1, maxdiff);
    end
    
catch ME
    fprintf('\n=== ERROR ===\n');
    fprintf('Error: %s\n', ME.message);
    rethrow(ME);
end

