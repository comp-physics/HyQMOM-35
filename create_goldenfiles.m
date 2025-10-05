% Script to create golden files for MPI regression testing

% Usage:
%   create_goldenfiles             % CI mode: creates 1 and 2 rank files
%   create_goldenfiles('local')    % Local mode: creates 4 and 8 rank files
%   create_goldenfiles('all')      % All: creates 1, 2, 4, and 8 rank files

% The CI-compatible files (1-2 ranks) are smaller and faster.
% Local files (4-8 ranks) provide comprehensive testing but require more resources.

function create_goldenfiles(mode)
    if nargin < 1
        mode = 'ci';  % Default to CI mode
    end
    
    fprintf('\n');
    
    switch lower(mode)
        case 'ci'
            fprintf('MPI Golden File Creation (CI Mode)');
            fprintf('Creates files for 1 and 2 MPI ranks');
            RANK_COUNTS = [1, 2];
        case 'local'
            fprintf('MPI Golden File Creation (Local Mode)');
            fprintf('Creates files for 4 and 8 MPI ranks');
            RANK_COUNTS = [4, 8];
        case 'all'
            fprintf('MPI Golden File Creation (All Modes)');
            fprintf('Creates files for 1, 2, 4, and 8 MPI ranks');
            RANK_COUNTS = [1, 2, 4, 8];
        otherwise
            error('Invalid mode. Use ''ci'', ''local'', or ''all''');
    end
    
    fprintf('\n');
    
    % Add src directory to path
    setup_paths();
    
    % Check required files
    required_files = {'main.m', 'src/setup_mpi_cartesian_2d.m', 'src/halo_exchange_2d.m'};
    for i = 1:length(required_files)
        if ~exist(required_files{i}, 'file')
            error('%s not found', required_files{i});
        end
    end
    
    % Check for Parallel Computing Toolbox
    has_pct = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    if ~has_pct
        error('Parallel Computing Toolbox is required for MPI golden file creation');
    end
    
    % Create goldenfiles directory if needed
    goldenfiles_dir = 'goldenfiles';
    if ~exist(goldenfiles_dir, 'dir')
        mkdir(goldenfiles_dir);
        fprintf('Created directory: %s\n\n', goldenfiles_dir);
    end
    
    % Parameters
    GOLDEN_TMAX = 0.1;
    
    % Determine grid sizes based on rank count
    % CI tests (1-2 ranks): 20x20 grid
    % Local tests (4-8 ranks): 40x40 grid
    NP_VALUES = zeros(size(RANK_COUNTS));
    for i = 1:length(RANK_COUNTS)
        if RANK_COUNTS(i) <= 2
            NP_VALUES(i) = 20;  % CI: smaller, faster
        else
            NP_VALUES(i) = 40;  % Local: larger, more comprehensive
        end
    end
    
    fprintf('Parameters:\n');
    fprintf('  tmax: %.2f\n', GOLDEN_TMAX);
    fprintf('  Rank counts: [%s]\n\n', num2str(RANK_COUNTS));
    
    % Verify grid sizes are valid
    fprintf('Grid configuration:\n');
    for i = 1:length(RANK_COUNTS)
        num_ranks = RANK_COUNTS(i);
        Np = NP_VALUES(i);
        [Px, Py] = mpi_utils('choose_grid', num_ranks);
        min_pts_x = floor(Np / Px);
        min_pts_y = floor(Np / Py);
        
        if min_pts_x < 10 || min_pts_y < 10
            error(['Grid too small for %d ranks (process grid %dx%d gives %d x %d pts/rank).\n' ...
                   'Need at least 10x10 per rank.'], ...
                   num_ranks, Px, Py, min_pts_x, min_pts_y);
        end
        
        fprintf('  %d rank(s): %dx%d grid, process grid %d x %d → %d x %d pts/rank\n', ...
                num_ranks, Np, Np, Px, Py, min_pts_x, min_pts_y);
    end
    fprintf('\n');
    
    %% Generate golden files
    all_success = true;
    for i = 1:length(RANK_COUNTS)
        num_ranks = RANK_COUNTS(i);
        Np = NP_VALUES(i);
        
        fprintf('  Generating golden file for %d rank(s) (Np=%d)\n', num_ranks, Np);
        
        try
            % Run simulation
            tic;
            fprintf('  Running simulation...\n');
            results = main(Np, GOLDEN_TMAX, false, num_ranks);
            elapsed = toc;
            
            fprintf('Simulation complete in %.1f seconds\n', elapsed);
            fprintf('Final time: %.6f\n', results.parameters.final_time);
            fprintf('Time steps: %d\n', results.parameters.time_steps);
            
            % Extract results (handle cell array from parallel pool)
            if iscell(results.moments.M)
                M_final = results.moments.M{1};
            else
                M_final = results.moments.M;
            end
            
            if iscell(results.grid)
                grid_final = results.grid{1};
            else
                grid_final = results.grid;
            end
            
            % Package golden data
            golden_data = struct();
            golden_data.M = M_final;
            golden_data.Np = Np;
            golden_data.tmax = GOLDEN_TMAX;
            golden_data.final_time = results.parameters.final_time;
            golden_data.time_steps = results.parameters.time_steps;
            golden_data.num_ranks = num_ranks;
            golden_data.xm = grid_final.xm;
            golden_data.ym = grid_final.ym;
            golden_data.creation_date = char(datetime("now", 'Format', 'dd-MMM-yyyy HH:mm:ss'));
            golden_data.matlab_version = version;
            
            % Get git branch if available
            [status, git_branch] = system('git rev-parse --abbrev-ref HEAD 2>/dev/null');
            if status == 0
                golden_data.git_branch = strtrim(git_branch);
            else
                golden_data.git_branch = 'unknown';
            end
            
            % Save golden file
            golden_filename = sprintf('goldenfile_mpi_%dranks_Np%d_tmax%d.mat', ...
                                      num_ranks, Np, round(GOLDEN_TMAX*1000));
            golden_path = fullfile(goldenfiles_dir, golden_filename);
            
            save(golden_path, 'golden_data', '-v7.3');
            
            file_info = dir(golden_path);
            fprintf('Saved: %s (%.1f KB)\n', golden_filename, file_info.bytes/1024);
            fprintf('Grid: %dx%d, %d ranks, %.3f final time\n\n', ...
                    Np, Np, num_ranks, results.parameters.final_time);
            
        catch ME
            fprintf('  ✗ ERROR: %s\n', ME.message);
            if ~isempty(ME.stack)
                fprintf('    at %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
            end
            all_success = false;
            fprintf('\n');
        end
    end
    
    %% Summary
    if all_success
        fprintf('GOLDEN FILE CREATION COMPLETE\n');
    else
        fprintf('GOLDEN FILE CREATION COMPLETED WITH ERRORS\n');
    end
    
    fprintf('Generated %d golden file(s) in %s/:\n', length(RANK_COUNTS), goldenfiles_dir);
    total_size = 0;
    for i = 1:length(RANK_COUNTS)
        num_ranks = RANK_COUNTS(i);
        Np = NP_VALUES(i);
        golden_filename = sprintf('goldenfile_mpi_%dranks_Np%d_tmax%d.mat', ...
                                  num_ranks, Np, round(GOLDEN_TMAX*1000));
        golden_path = fullfile(goldenfiles_dir, golden_filename);
        if exist(golden_path, 'file')
            finfo = dir(golden_path);
            total_size = total_size + finfo.bytes;
            fprintf('✓ %s (%.1f KB)\n', golden_filename, finfo.bytes/1024);
        else
            fprintf('✗ %s (FAILED)\n', golden_filename);
        end
    end
    
    fprintf('\nTotal size: %.1f KB\n\n', total_size/1024);
    
    % Mode-specific next steps
    fprintf('Next steps:\n');
    switch lower(mode)
        case 'ci'
            fprintf('  1. Run: cd tests && runtests(''test_mpi_goldenfile'')\n');
            fprintf('  2. These files are used by CI (1-2 ranks only)\n');
            fprintf('  3. Commit these files to git if intentionally changed\n');
        case 'local'
            fprintf('  1. Run: cd tests && runtests(''test_mpi_local_only'')\n');
            fprintf('  2. These files are for local testing (4-8 ranks)\n');
            fprintf('  3. Commit these files to git:\n');
            fprintf('     git add goldenfiles/goldenfile_mpi_4ranks*.mat\n');
            fprintf('     git add goldenfiles/goldenfile_mpi_8ranks*.mat\n');
            fprintf('     git commit -m "Update local MPI golden files"\n');
        case 'all'
            fprintf('  1. Run: cd tests && runtests  % Runs all tests\n');
            fprintf('  2. Commit updated files to git if intentionally changed\n');
    end
    
    fprintf('\nNote: Golden files are committed to version control.\n');
    fprintf('      CI tests auto-skip local-only tests (4+ ranks).\n');
    fprintf('\n');
end
