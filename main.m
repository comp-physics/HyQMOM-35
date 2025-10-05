function [results] = main(varargin)
% Main solver for 3D HyQMOM with MPI-parallel domain decomposition
% Parameters:
%   Np           - GLOBAL grid size (total points in each direction)
%   tmax         - Final simulation time
%   enable_plots - Enable/disable plotting (default: false)
%   num_workers  - Number of MPI ranks/workers (default: 2)
%   enable_profile - Enable MPI profiling (default: false)
% Usage:
%   main()                           % Run with defaults
%   main(Np, tmax)                   % Override Np and tmax
%   main(Np, tmax, enable_plots)     % Override plotting
%   main(Np, tmax, enable_plots, num_workers) % Specify number of MPI ranks
%   main(Np, tmax, enable_plots, num_workers, enable_profile) % Enable profiling
% Examples:
%   main()                    % Default: Np=20 (global), tmax=0.05, 2 workers
%   main(40, 0.1, false, 2)   % 40×40 GLOBAL grid, 2 MPI ranks (each gets 40×20)
%   main(40, 0.1, false, 4)   % 40×40 GLOBAL grid, 4 MPI ranks (each gets 20×20)
%   main(40, 0.1, false, 4, true) % Same as above with MPI profiling enabled
% Note: Np is the TOTAL grid size. It will be decomposed into subdomains.
%       Each rank must have at least 10×10 interior points.

% Add src directory to path
script_dir = fileparts(mfilename('fullpath'));
setup_paths(script_dir);

% Parse input arguments
defaults = struct('Np', 140, 'tmax', 0.02, 'enable_plots', false, 'num_workers', 6, 'enable_profile', false);
if nargin == 0
    Np = defaults.Np;
    tmax = defaults.tmax;
    enable_plots = defaults.enable_plots;
    num_workers = defaults.num_workers;
    enable_profile = defaults.enable_profile;
elseif nargin == 1
    Np = varargin{1};
    tmax = defaults.tmax;
    enable_plots = defaults.enable_plots;
    num_workers = defaults.num_workers;
    enable_profile = defaults.enable_profile;
elseif nargin == 2
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = defaults.enable_plots;
    num_workers = defaults.num_workers;
    enable_profile = defaults.enable_profile;
elseif nargin == 3
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = varargin{3};
    num_workers = defaults.num_workers;
    enable_profile = defaults.enable_profile;
elseif nargin == 4
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = varargin{3};
    num_workers = varargin{4};
    enable_profile = defaults.enable_profile;
else
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = varargin{3};
    num_workers = varargin{4};
    enable_profile = varargin{5};
end

% Physical parameters
Kn = 1.0;
Ma = 0.0;
flag2D = 0;

% Spatial and temporal discretization
CFL = 0.5;
dx = 1.0 / Np;
dy = 1.0 / Np;

N = 4;
Nmom = 35;
Nmom5 = 21;
nnmax = 2e7;
dtmax = Kn;

% Correlation coefficients
r110 = 0.0;
r101 = 0.0;
r011 = 0.0;

% Initial conditions
T = 1.0;
rhol = 1.0;
rhor = 0.01;

% Validate grid size for MPI decomposition
% Requirement: minimum ~10 points per rank in each direction
[Px, Py] = mpi_utils('choose_grid', num_workers);
min_points_x = floor(Np / Px);
min_points_y = floor(Np / Py);
min_points = min(min_points_x, min_points_y);

if min_points < 10
    error(['Grid too small for %d workers. With Np=%d, process grid %dx%d gives ' ...
           'only %d points/rank (minimum 10 required). Use fewer workers or larger Np.'], ...
           num_workers, Np, Px, Py, min_points);
end

% Start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', num_workers);
elseif pool.NumWorkers ~= num_workers
    delete(pool);
    parpool('local', num_workers);
end

% Initialize MPI profiler if enabled
if enable_profile
    fprintf('MPI profiling enabled. Starting mpiprofile...\n');
    mpiprofile on;
end

% Pass script directory to workers (determine before spmd)
script_dir_for_workers = script_dir;

% MPI parallel execution
spmd
    % Workers need src/ directory in their path
    setup_paths(script_dir_for_workers);
    
    % Define all constants directly inside spmd block
    halo = 2;  % Required for pas_HLL stencil
    bc = struct('type', 'copy');
    
    % Physical parameters (must match values outside spmd)
    Kn_worker = 1.0;
    Ma_worker = 0.0;
    flag2D_worker = 0;
    CFL_worker = 0.5;
    dx_worker = 1.0 / Np;
    dy_worker = 1.0 / Np;
    Nmom_worker = 35;
    Nmom5_worker = 21;
    nnmax_worker = 2e7;
    dtmax_worker = 1.0;  % Same as Kn
    r110_worker = 0.0;
    r101_worker = 0.0;
    r011_worker = 0.0;
    T_worker = 1.0;
    rhol_worker = 1.0;
    rhor_worker = 0.01;
    
    % Grid setup (each worker needs this)
    grid = grid_utils('setup', Np, -0.5, 0.5, -0.5, 0.5);
    M_global = grid_utils('crossing_jets_ic', Np, Nmom_worker, rhol_worker, rhor_worker, Ma_worker, T_worker, r110_worker, r101_worker, r011_worker);
    
    % Setup domain decomposition
    decomp = setup_mpi_cartesian_2d(Np, halo);
    nx = decomp.local_size(1);
    ny = decomp.local_size(2);
    
    % Allocate local arrays with halos
    M = zeros(nx+2*halo, ny+2*halo, Nmom_worker);
    Mnp = M;
    Fx = zeros(nx+2*halo, ny+2*halo, Nmom_worker);
    Fy = zeros(nx+2*halo, ny+2*halo, Nmom_worker);
    
    % Wave speed arrays (interior only, no halos needed)
    vpxmin = zeros(nx, ny);
    vpxmax = zeros(nx, ny);
    vpymin = zeros(nx, ny);
    vpymax = zeros(nx, ny);
    v5xmin = zeros(nx, ny);
    v5xmax = zeros(nx, ny);
    v5ymin = zeros(nx, ny);
    v5ymax = zeros(nx, ny);
    v6xmin = zeros(nx, ny);
    v6xmax = zeros(nx, ny);
    v6ymin = zeros(nx, ny);
    v6ymax = zeros(nx, ny);
    
    % Scatter initial conditions to local subdomains
    i0i1 = decomp.istart_iend;
    j0j1 = decomp.jstart_jend;
    M(halo+1:halo+nx, halo+1:halo+ny, :) = M_global(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :);
    
    % Initial halo exchange
    M = halo_exchange_2d(M, decomp, bc);
    
    % Time evolution
    t = 0.0;
    nn = 0;
    
    while t < tmax && nn < nnmax_worker
        nn = nn + 1;
        step_start_time = tic;  % Start timing this timestep
        
        % Compute fluxes and wave speeds for interior cells
        for i = 1:nx
            for j = 1:ny
                % Access with halo offset
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D_worker, Ma_worker);
                [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D_worker, Ma_worker);
                [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D_worker, Ma_worker);
                
                Fx(ih, jh, :) = Mx;
                Fy(ih, jh, :) = My;
                Mnp(ih, jh, :) = Mr;
                
                [~, v5xmin(i,j), v5xmax(i,j)] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
                [~, v5ymin(i,j), v5ymax(i,j)] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
                
                vpxmin(i,j) = min(v5xmin(i,j), v6xmin(i,j));
                vpxmax(i,j) = max(v5xmax(i,j), v6xmax(i,j));
                vpymin(i,j) = min(v5ymin(i,j), v6ymin(i,j));
                vpymax(i,j) = max(v5ymax(i,j), v6ymax(i,j));
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Exchange M, Fx, Fy from neighbors
        M = halo_exchange_2d(M, decomp, bc);
        Fx = halo_exchange_2d(Fx, decomp, bc);
        Fy = halo_exchange_2d(Fy, decomp, bc);
        
        % Compute fluxes and wave speeds in halo cells for pas_HLL stencil
        [Fx, Fy, vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext] = ...
            compute_halo_fluxes_and_wavespeeds(M, Fx, Fy, vpxmin, vpxmax, vpymin, vpymax, ...
                                               nx, ny, halo, flag2D_worker, Ma_worker);
        
        % Global reduction for time step (all ranks need same dt)
        vmax_local = max([abs(vpxmax(:)); abs(vpxmin(:)); abs(vpymax(:)); abs(vpymin(:))]);
        vmax = max(gcat(vmax_local, 1));  % Global max across all ranks
        dt = min(CFL_worker*dx_worker/vmax, dtmax_worker);
        dt = min(dt, tmax-t);
        t = t + dt;
        
        % X-direction flux update with processor boundary handling
        Mnpx = apply_flux_update(M, Fx, vpxmin, vpxmax, vpxmin_ext, vpxmax_ext, ...
                                  nx, ny, halo, dt, dx_worker, decomp, 1);
        
        % Y-direction flux update with processor boundary handling
        Mnpy = apply_flux_update(M, Fy, vpymin, vpymax, vpymin_ext, vpymax_ext, ...
                                  nx, ny, halo, dt, dy_worker, decomp, 2);
        
        % Combine updates (Strang splitting) - INTERIOR ONLY
        % Mnpx and Mnpy halos contain stale M values, only interior was updated by pas_HLL
        M(halo+1:halo+nx, halo+1:halo+ny, :) = ...
            Mnpx(halo+1:halo+nx, halo+1:halo+ny, :) + ...
            Mnpy(halo+1:halo+nx, halo+1:halo+ny, :) - ...
            M(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Exchange halos before realizability enforcement
        M = halo_exchange_2d(M, decomp, bc);
        
        % Enforce realizability
        for i = 1:nx
            for j = 1:ny
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D_worker, Ma_worker);
                [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D_worker, Ma_worker);
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D_worker, Ma_worker);
                Mnp(ih, jh, :) = Mr;
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Apply BGK collision
        for i = 1:nx
            for j = 1:ny
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                MMC = collision35(MOM, dt, Kn_worker);
                Mnp(ih, jh, :) = MMC;
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Exchange halos for next iteration
        M = halo_exchange_2d(M, decomp, bc);
        
        % Compute MaxDiff for symmetry check over GLOBAL domain
        % Gather interior data from all ranks to rank 1
        M_interior_local = M(halo+1:halo+nx, halo+1:halo+ny, :);
        
        if spmdIndex == 1
            % Rank 1: collect from all workers to reconstruct global array
            M_global_temp = zeros(Np, Np, Nmom);
            M_global_temp(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :) = M_interior_local;
            
            % Receive from all other ranks
            for src = 2:spmdSize
                data_packet = spmdReceive(src);
                blk = data_packet{1};
                i_range = data_packet{2};
                j_range = data_packet{3};
                M_global_temp(i_range(1):i_range(2), j_range(1):j_range(2), :) = blk;
            end
            
            % Now compute MaxDiff on full global domain
            [~, MaxDiff] = diagnostics('test_symmetry', M_global_temp, Np);
        else
            % All other ranks: send to rank 1
            data_packet = {M_interior_local, i0i1, j0j1};
            spmdSend(data_packet, 1);
            MaxDiff = [];  % Will be broadcast from rank 1
        end
        
        % Compute timing on all ranks
        step_time = toc(step_start_time);
        
        % Compute time per grid point for this rank
        local_grid_points = nx * ny;
        time_per_point_local = step_time / local_grid_points;
        
        % Gather max time per point across all ranks
        max_time_per_point = gop(@max, time_per_point_local);
        
        % Print timestep timing and MaxDiff (only from rank 1)
        if spmdIndex == 1
            fprintf('Step %4d: t = %.6f, dt = %.6e, max s/pt = %.6e s, MaxDiff = %.3e\n', ...
                    nn, t, dt, max_time_per_point, max(abs(MaxDiff)));
        end
    end
    
    % Gather results to rank 1 using asynchronous send/receive
    M_interior = M(halo+1:halo+nx, halo+1:halo+ny, :);
    
    if spmdIndex == 1
        % Rank 1: gather from all ranks using mpi_utils
        M_final = mpi_utils('gather_M', M_interior, i0i1, j0j1, Np, Nmom);
        final_time = t;
        time_steps = nn;
    else
        % All other ranks: send to rank 1 using mpi_utils
        mpi_utils('send_M', M_interior, i0i1, j0j1, 1);
        M_final = [];
        final_time = [];
        time_steps = [];
    end
end

% Collect and display MPI profiling results if enabled
if enable_profile
    fprintf('\nCollecting MPI profile data...\n');
    profile_stats = mpiprofile('info');
    
    % Save profile data to file
    profile_filename = sprintf('mpi_profile_%dranks_Np%d.mat', num_workers, Np);
    save(profile_filename, 'profile_stats');
    fprintf('Profile data saved to: %s\n', profile_filename);
    fprintf('To view later, use: load(''%s''); mpiprofile(''viewer'', profile_stats)\n', profile_filename);
    
    % Open profile viewer
    fprintf('Opening MPI profile viewer...\n');
    mpiprofile viewer;
    
    % Display summary statistics
    fprintf('\n=== MPI Profile Summary ===\n');
    for w = 1:length(profile_stats)
        if ~isempty(profile_stats(w).FunctionTable)
            worker_time = sum([profile_stats(w).FunctionTable.TotalTime]);
            fprintf('Worker %d: Total time = %.2f s\n', w, worker_time);
        end
    end
    fprintf('===========================\n\n');
end

% Extract from composite
if isa(M_final, 'Composite')
    M_final = M_final{1};
    final_time = final_time{1};
    time_steps = time_steps{1};
    grid_worker = grid{1};  % Extract grid from first worker
else
    grid_worker = grid;
end

% Compute derived quantities on gathered result
[C, S] = grid_moment_processor(M_final, @M2CS4_35);
[M5, C5, S5] = grid_moment_processor(M_final, @Moments5_3D);

% Plot results if requested using comprehensive plotting functions
if enable_plots
    eig_data = grid_eigenvalues(M_final, Np, Nmom);
    simulation_plots('final', grid_worker.xm, grid_worker.ym, M_final, C, S, M5, C5, S5, ...
                     Np, eig_data, enable_plots);
    drawnow;
end

% Build results structure
if nargout > 0
    results = struct();
    results.parameters.Np = Np;
    results.parameters.tmax = tmax;
    results.parameters.enable_plots = enable_plots;
    results.parameters.num_workers = num_workers;
    results.parameters.enable_profile = enable_profile;
    results.parameters.Kn = Kn;
    results.parameters.Ma = Ma;
    results.parameters.CFL = CFL;
    results.parameters.Nmom = Nmom;
    results.parameters.N = N;
    results.parameters.final_time = final_time;
    results.parameters.time_steps = time_steps;
    
    results.grid.x = grid_worker.x;
    results.grid.y = grid_worker.y;
    results.grid.xm = grid_worker.xm;
    results.grid.ym = grid_worker.ym;
    results.grid.dx = grid_worker.dx;
    results.grid.dy = grid_worker.dy;
    
    results.moments.M = M_final;
    results.moments.C = C;
    results.moments.S = S;
    results.moments.M5 = M5;
    results.moments.C5 = C5;
    results.moments.S5 = S5;
    
    results.filename = sprintf('mpi_%dranks_Np%d_tmax%g.mat', num_workers, Np, tmax);
end

end

