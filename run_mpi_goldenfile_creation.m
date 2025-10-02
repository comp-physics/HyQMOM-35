% Script to create golden files for MPI regression testing
% Creates golden files for 1, 2, and 4 MPI ranks

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║                                                              ║\n');
fprintf('║     MPI Golden File Creation                                ║\n');
fprintf('║                                                              ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

% Add src directory to path
addpath('src');

% Check required files
required_files = {'main.m', 'src/setup_mpi_cartesian_2d.m', 'src/halo_exchange_2d.m'};
for i = 1:length(required_files)
    if ~exist(required_files{i}, 'file')
        error('%s not found', required_files{i});
    end
end

% Create goldenfiles directory if needed
goldenfiles_dir = 'goldenfiles';
if ~exist(goldenfiles_dir, 'dir')
    mkdir(goldenfiles_dir);
    fprintf('Created directory: %s\n\n', goldenfiles_dir);
end

% Parameters (Np=20 required for halo=2 with 4 ranks: 20/2=10 pts/rank minimum)
GOLDEN_NP = 20;
GOLDEN_TMAX = 0.1;
RANK_COUNTS = [1, 2, 4];

fprintf('Golden file parameters:\n');
fprintf('  Grid size: %d × %d\n', GOLDEN_NP, GOLDEN_NP);
fprintf('  Final time: %.3f\n', GOLDEN_TMAX);
fprintf('  Rank counts: %s\n\n', mat2str(RANK_COUNTS));

all_success = true;

for i = 1:length(RANK_COUNTS)
    num_ranks = RANK_COUNTS(i);
    
    fprintf('═══════════════════════════════════════════════════════════\n');
    fprintf('Creating golden file for %d rank(s)...\n', num_ranks);
    fprintf('═══════════════════════════════════════════════════════════\n\n');
    
    try
        % Run MPI simulation using unified main() interface
        fprintf('Running MPI simulation with %d rank(s)...\n', num_ranks);
        tic;
        results = main(GOLDEN_NP, GOLDEN_TMAX, false, false, true, num_ranks);
        elapsed_time = toc;
        
        fprintf('Simulation completed in %.2f seconds\n', elapsed_time);
        fprintf('  Final time: %.6f\n', results.parameters.final_time);
        fprintf('  Time steps: %d\n\n', results.parameters.time_steps);
        
        % Create golden data structure
        golden_data = results;
        golden_data.parameters.elapsed_time = elapsed_time;
        
        % Add metadata
        golden_data.metadata.creation_date = datestr(now);
        golden_data.metadata.matlab_version = version;
        golden_data.metadata.description = sprintf('MPI golden file: Np=%d, tmax=%.3f, ranks=%d', ...
                                                    GOLDEN_NP, GOLDEN_TMAX, num_ranks);
        golden_data.metadata.num_ranks = num_ranks;
        golden_data.metadata.git_branch = 'mpi';
        
        % Save golden file
        golden_filename = fullfile(goldenfiles_dir, ...
                                  sprintf('goldenfile_mpi_%dranks_Np%d_tmax%03d.mat', ...
                                          num_ranks, GOLDEN_NP, round(GOLDEN_TMAX*1000)));
        save(golden_filename, 'golden_data');
        
        fprintf('✓ Golden file saved: %s\n', golden_filename);
        fprintf('  File size: %.2f KB\n', dir(golden_filename).bytes / 1024);
        
        % Verify symmetry
        fprintf('\nVerifying symmetry:\n');
        M = golden_data.moments.M;
        max_asym = 0;
        for mom = 1:5
            asym = max(abs(M(:,:,mom) - flipud(M(:,:,mom))), [], 'all');
            fprintf('  M%d asymmetry: %.3e\n', mom-1, asym);
            max_asym = max(max_asym, asym);
        end
        
        if max_asym < 1e-12
            fprintf('  ✓ Symmetry preserved\n');
        else
            fprintf('  ⚠ Warning: symmetry degraded (max: %.3e)\n', max_asym);
        end
        
        fprintf('\n');
        
    catch ME
        fprintf('\n✗ ERROR creating golden file for %d rank(s):\n', num_ranks);
        fprintf('  %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  Stack trace:\n');
            for j = 1:min(3, length(ME.stack))
                fprintf('    %s (line %d)\n', ME.stack(j).name, ME.stack(j).line);
            end
        end
        all_success = false;
        fprintf('\n');
    end
end

% Summary
fprintf('═══════════════════════════════════════════════════════════\n');
if all_success
    fprintf('✓ SUCCESS: All MPI golden files created\n');
    fprintf('\nCreated files:\n');
    for i = 1:length(RANK_COUNTS)
        filename = sprintf('goldenfile_mpi_%dranks_Np%d_tmax%03d.mat', ...
                          RANK_COUNTS(i), GOLDEN_NP, round(GOLDEN_TMAX*1000));
        filepath = fullfile(goldenfiles_dir, filename);
        if exist(filepath, 'file')
            fprintf('  ✓ %s (%.2f KB)\n', filename, dir(filepath).bytes / 1024);
        end
    end
    
    fprintf('\nNext steps:\n');
    fprintf('  1. Run: runtests(''tests/test_mpi_goldenfile.m'')\n');
    fprintf('  2. Compare MPI results with serial golden file\n');
    fprintf('  3. Verify consistency across different rank counts\n');
else
    fprintf('⚠ WARNING: Some golden files failed to create\n');
    fprintf('Check error messages above for details\n');
end
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('\n');

