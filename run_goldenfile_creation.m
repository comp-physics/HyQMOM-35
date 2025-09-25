% Script to create golden files for regression testing
% This script runs the simulation with fixed parameters and saves the results

fprintf('=== GOLDEN FILE CREATION WRAPPER ===\n');
fprintf('This script will create a golden file by running the simulation with:\n');
fprintf('  - Np = 6 (grid points)\n');
fprintf('  - tmax = 0.02 (final time)\n');
fprintf('  - enable_plots = false (no plotting)\n\n');

% Check if we're in the right directory
required_files = {'main_2Dcrossing_3DHyQMOM35.m', 'simulation_plots.m'};
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

% Check if golden file already exists
golden_filename = fullfile(goldenfiles_dir, 'goldenfile_Np6_tmax020.mat');
if exist(golden_filename, 'file')
    response = input(sprintf('Golden file %s already exists. Overwrite? (y/n): ', golden_filename), 's');
    if ~strcmpi(response, 'y')
        fprintf('Aborted by user.\n');
        return;
    end
    fprintf('Overwriting existing golden file...\n');
end

% Define golden file parameters (these will override the main script values)
GOLDEN_NP = 6;
GOLDEN_TMAX = 0.02;
GOLDEN_ENABLE_PLOTS = false;

% Run the golden file creation
try
    fprintf('Golden file parameters: Np = %d, tmax = %.3f, enable_plots = %s\n', ...
        GOLDEN_NP, GOLDEN_TMAX, mat2str(GOLDEN_ENABLE_PLOTS));
    fprintf('Running simulation with enforced parameters...\n');
    tic;
    
    % Execute the main simulation with parameter overrides
    results = main_2Dcrossing_3DHyQMOM35(GOLDEN_NP, GOLDEN_TMAX, GOLDEN_ENABLE_PLOTS);
    
    elapsed_time = toc;
    fprintf('Simulation completed in %.2f seconds\n', elapsed_time);
    
    % Create golden file structure from returned results
    golden_data = results;
    
    % Add elapsed time to parameters
    golden_data.parameters.elapsed_time = elapsed_time;
    
    % Add metadata
    golden_data.metadata.creation_date = datestr(now);
    golden_data.metadata.matlab_version = version;
    golden_data.metadata.description = sprintf('Golden file for main_2Dcrossing_3DHyQMOM35.m with Np=%d, tmax=%.3f', GOLDEN_NP, GOLDEN_TMAX);
    
    % Save to file in goldenfiles directory
    golden_filename = fullfile(goldenfiles_dir, sprintf('goldenfile_Np%d_tmax%03d.mat', GOLDEN_NP, round(GOLDEN_TMAX*1000)));
    save(golden_filename, 'golden_data');
    
    fprintf('Golden file saved as: %s\n', golden_filename);
    fprintf('File size: %.2f KB\n', dir(golden_filename).bytes / 1024);
    
    % Display summary of saved data
    fprintf('\nSaved data summary:\n');
    fprintf('  - Parameters: %d fields\n', length(fieldnames(golden_data.parameters)));
    fprintf('  - Grid points: %dx%d\n', length(golden_data.grid.xm), length(golden_data.grid.ym));
    fprintf('  - Moment arrays: %d fields\n', length(fieldnames(golden_data.moments)));
    if isfield(golden_data, 'eigenvalues')
        fprintf('  - Eigenvalue data: %d fields\n', length(fieldnames(golden_data.eigenvalues)));
    end
    if isfield(golden_data, 'velocities')
        fprintf('  - Velocity bounds: %d fields\n', length(fieldnames(golden_data.velocities)));
    end
    
    fprintf('\n=== SUCCESS ===\n');
    fprintf('Golden file created successfully!\n');
    fprintf('You can now use validate_against_goldenfile.m to test future changes.\n');
    
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
