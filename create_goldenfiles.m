% Script to create golden files for MPI regression testing
% Creates golden files for 1, 2, 3, and 4 MPI ranks
% Each configuration uses 20 grid points per rank per dimension

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║                                                              ║\n');
fprintf('║     MPI Golden File Creation (CI-Compatible)                ║\n');
fprintf('║     Creates golden files for 1 and 2 MPI ranks              ║\n');
fprintf('║     (CI environment limited to 2 workers)                   ║\n');
fprintf('║                                                              ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

% Add src directory to path
addpath('src');
addpath('src/autogen');

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

% Parameters: 20 grid points per rank per dimension
POINTS_PER_RANK = 20;
GOLDEN_TMAX = 0.1;
% Create golden files for 1 and 2 ranks (CI limitation: max 2 workers)
RANK_COUNTS = [1, 2];

% Calculate Np for each rank count to ensure >= 10 pts/rank in BOTH directions
% Np must be divisible by both Px and Py, with at least 10 points per rank
NP_VALUES = zeros(size(RANK_COUNTS));
for i = 1:length(RANK_COUNTS)
    [Px, Py] = mpi_utils('choose_grid', RANK_COUNTS(i));
    % Calculate Np to give at least POINTS_PER_RANK in each direction
    % Np = max(Px, Py) * POINTS_PER_RANK ensures all ranks get POINTS_PER_RANK
    Np_candidate = max(Px, Py) * POINTS_PER_RANK;
    % But we need Np / Px >= 10 AND Np / Py >= 10
    min_for_x = Px * 10;  % Minimum to give 10 pts/rank in x
    min_for_y = Py * 10;  % Minimum to give 10 pts/rank in y
    Np_min = max(min_for_x, min_for_y);
    NP_VALUES(i) = max(Np_candidate, Np_min);
end

fprintf('Golden file parameters:\n');
fprintf('  Grid points per rank: %d × %d\n', POINTS_PER_RANK, POINTS_PER_RANK);
fprintf('  Final time: %.3f\n', GOLDEN_TMAX);
fprintf('  Rank counts: %s\n\n', mat2str(RANK_COUNTS));

% Show grid sizes
fprintf('Grid sizes (to achieve ~20 pts/rank minimum in each direction):\n');
for i = 1:length(RANK_COUNTS)
    num_ranks = RANK_COUNTS(i);
    [Px, Py] = choose_process_grid(num_ranks);
    Np = NP_VALUES(i);
    pts_per_rank_x = Np / Px;
    pts_per_rank_y = Np / Py;
    fprintf('  %d rank(s): %dx%d process grid → %d×%d grid (%.1f×%.1f pts/rank)\n', ...
            num_ranks, Px, Py, Np, Np, pts_per_rank_x, pts_per_rank_y);
end
fprintf('\n');

all_success = true;

for i = 1:length(RANK_COUNTS)
    num_ranks = RANK_COUNTS(i);
    
    % Get pre-calculated grid size for this rank count
    Np = NP_VALUES(i);
    [Px, Py] = choose_process_grid(num_ranks);
    
    fprintf('═══════════════════════════════════════════════════════════\n');
    fprintf('Creating golden file for %d rank(s) (%d×%d grid)...\n', num_ranks, Np, Np);
    fprintf('═══════════════════════════════════════════════════════════\n\n');
    
    try
        % Run MPI simulation
        fprintf('Running MPI simulation with %d rank(s)...\n', num_ranks);
        tic;
        results = main(Np, GOLDEN_TMAX, false, num_ranks);
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
        golden_data.metadata.description = sprintf('MPI golden file: %d ranks, Np=%d (%d pts/rank), tmax=%.3f', ...
                                                    num_ranks, Np, POINTS_PER_RANK, GOLDEN_TMAX);
        golden_data.metadata.num_ranks = num_ranks;
        golden_data.metadata.points_per_rank = POINTS_PER_RANK;
        golden_data.metadata.git_branch = 'mpi';
        
        % Save golden file
        golden_filename = fullfile(goldenfiles_dir, ...
                                  sprintf('goldenfile_mpi_%dranks_Np%d_tmax%03d.mat', ...
                                          num_ranks, Np, round(GOLDEN_TMAX*1000)));
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
        num_ranks = RANK_COUNTS(i);
        Np = NP_VALUES(i);
        filename = sprintf('goldenfile_mpi_%dranks_Np%d_tmax%03d.mat', ...
                          num_ranks, Np, round(GOLDEN_TMAX*1000));
        filepath = fullfile(goldenfiles_dir, filename);
        if exist(filepath, 'file')
            fprintf('  ✓ %s (%.2f KB)\n', filename, dir(filepath).bytes / 1024);
        end
    end
    
    fprintf('\nNext steps:\n');
    fprintf('  1. Run: runtests(''tests/test_mpi_goldenfile.m'')\n');
    fprintf('  2. Tests are CI-compatible (max 2 workers)\n');
else
    fprintf('⚠ WARNING: Some golden files failed to create\n');
    fprintf('Check error messages above for details\n');
end
fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('\n');

