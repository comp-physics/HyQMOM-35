function [results] = main(varargin)
% Main solver for 3D HyQMOM with MPI-parallel domain decomposition
% Parameters:
%   Np           - GLOBAL grid size (total points in x and y directions)
%   tmax         - Final simulation time
%   enable_plots - Enable/disable plotting (default: false)
%   num_workers  - Number of MPI ranks/workers (default: 6)
%   enable_profile - Enable MPI profiling (default: false)
%   symmetry_check_interval - Check symmetry every N steps (default: 1)
%                             Set to 100 for large problems to improve performance
%   enable_memory_tracking - Enable memory usage tracking (default: true)
%   Nz           - Grid size in z direction (default: 4, no MPI decomposition in z)
% Usage:
%   main()                           % Run with defaults
%   main(Np, tmax)                   % Override Np and tmax
%   main(Np, tmax, enable_plots)     % Override plotting
%   main(Np, tmax, enable_plots, num_workers) % Specify number of MPI ranks
%   main(Np, tmax, enable_plots, num_workers, enable_profile) % Enable profiling
%   main(Np, tmax, enable_plots, num_workers, enable_profile, symmetry_check_interval) % Set check interval
%   main(Np, tmax, enable_plots, num_workers, enable_profile, symmetry_check_interval, enable_memory_tracking) % Control memory tracking
%   main(Np, tmax, enable_plots, num_workers, enable_profile, symmetry_check_interval, enable_memory_tracking, Nz) % Set z grid size
% Examples:
%   main()                    % Default: Np=120 (global), Nz=4, tmax=0.02, 6 workers
%   main(40, 0.1, false, 2)   % 40x40x4 GLOBAL grid, 2 MPI ranks (each gets 40x20x4)
%   main(40, 0.1, false, 4)   % 40x40x4 GLOBAL grid, 4 MPI ranks (each gets 20x20x4)
%   main(40, 0.1, false, 4, true) % Same as above with MPI profiling enabled
%   main(200, 0.1, false, 4, false, 100) % Check symmetry every 100 steps (faster!)
%   main(40, 0.1, false, 4, false, 1, false, 8) % 40x40x8 grid
% Note: Np is the total grid size in x and y. It will be decomposed into subdomains.
%       Nz is NOT decomposed (all ranks have full z extent).
%       Each rank must have at least 10x10 interior points in x-y.

% Add src directory to path
script_dir = fileparts(mfilename('fullpath'));
setup_paths(script_dir);

% Parse input arguments with clean helper
defaults = struct('Np', 20, 'tmax', 0.1, 'enable_plots', false, ...
                  'num_workers', 1, 'enable_profile', false, ...
                  'symmetry_check_interval', 10, 'enable_memory_tracking', false, ...
                  'Nz', 4, 'homogeneous_z', true);
params = parse_main_args(varargin, defaults);

% Validate grid size for MPI decomposition
% Requirement: minimum ~10 points per rank in each direction
[Px, Py] = mpi_utils('choose_grid', params.num_workers);
min_points_x = floor(params.Np / Px);
min_points_y = floor(params.Np / Py);
min_points = min(min_points_x, min_points_y);

if min_points < 10
    error(['Grid too small for %d workers. With Np=%d, process grid %dx%d gives ' ...
           'only %d points/rank (minimum 10 required). Use fewer workers or larger Np.'], ...
           params.num_workers, params.Np, Px, Py, min_points);
end

% Start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', params.num_workers);
elseif pool.NumWorkers ~= params.num_workers
    delete(pool);
    pool = parpool('local', params.num_workers);
end

% Initialize MPI profiler if enabled
if params.enable_profile
    profiling_utils('start');
end

% Run parallel simulation (core SPMD loop)
[M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir);

% Collect and display MPI profiling results if enabled
if params.enable_profile
    profiling_utils('report', params.num_workers, params.Np);
end

% Compute derived quantities on gathered result
[C, S] = grid_moment_processor(M_final, @M2CS4_35);
[M5, C5, S5] = grid_moment_processor(M_final, @Moments5_3D);

% Plot results if requested using comprehensive plotting functions
if params.enable_plots
    eig_data = grid_eigenvalues(M_final, params.Np, params.Nmom);
    simulation_plots('final', grid_out.xm, grid_out.ym, M_final, C, S, M5, C5, S5, ...
                     params.Np, eig_data, params.enable_plots);
    drawnow;
end

% Build results structure
if nargout > 0
    results = struct();
    results.parameters.Np = params.Np;
    results.parameters.tmax = params.tmax;
    results.parameters.enable_plots = params.enable_plots;
    results.parameters.num_workers = params.num_workers;
    results.parameters.enable_profile = params.enable_profile;
    results.parameters.Kn = params.Kn;
    results.parameters.Ma = params.Ma;
    results.parameters.CFL = params.CFL;
    results.parameters.Nmom = params.Nmom;
    results.parameters.N = params.N;
    results.parameters.symmetry_check_interval = params.symmetry_check_interval;
    results.parameters.enable_memory_tracking = params.enable_memory_tracking;
    results.parameters.final_time = final_time;
    results.parameters.time_steps = time_steps;
    
    results.grid.x = grid_out.x;
    results.grid.y = grid_out.y;
    results.grid.z = grid_out.z;
    results.grid.xm = grid_out.xm;
    results.grid.ym = grid_out.ym;
    results.grid.zm = grid_out.zm;
    results.grid.dx = grid_out.dx;
    results.grid.dy = grid_out.dy;
    results.grid.dz = grid_out.dz;
    
    results.moments.M = M_final;
    results.moments.C = C;
    results.moments.S = S;
    results.moments.M5 = M5;
    results.moments.C5 = C5;
    results.moments.S5 = S5;
    
    results.filename = sprintf('mpi_%dranks_Np%d_tmax%g.mat', params.num_workers, params.Np, params.tmax);
end

end