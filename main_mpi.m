function [results] = main_mpi(varargin)
% MPI-parallel version of main solver using domain decomposition
% Same interface as main.m but runs with MPI parallelization
%
% Usage:
%   main_mpi()                           % Run with defaults
%   main_mpi(Np, tmax)                   % Override Np and tmax
%   main_mpi(Np, tmax, enable_plots)     % Override plotting
%   main_mpi(Np, tmax, enable_plots, num_workers) % Specify number of MPI ranks
%
% Examples:
%   main_mpi()                    % Default: Np=10, tmax=0.1, 4 workers
%   main_mpi(10, 0.1, false, 2)   % 10x10 grid, 2 MPI ranks

% Add src directory to path
script_dir = fileparts(mfilename('fullpath'));
src_dir = fullfile(script_dir, 'src');
if exist(src_dir, 'dir')
    addpath(src_dir);
end

% Parse input arguments
defaults = struct('Np', 10, 'tmax', 0.1, 'enable_plots', false, 'num_workers', 4);
if nargin == 0
    Np = defaults.Np;
    tmax = defaults.tmax;
    enable_plots = defaults.enable_plots;
    num_workers = defaults.num_workers;
elseif nargin == 1
    Np = varargin{1};
    tmax = defaults.tmax;
    enable_plots = defaults.enable_plots;
    num_workers = defaults.num_workers;
elseif nargin == 2
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = defaults.enable_plots;
    num_workers = defaults.num_workers;
elseif nargin == 3
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = varargin{3};
    num_workers = defaults.num_workers;
else
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = varargin{3};
    num_workers = varargin{4};
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

% Start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', num_workers);
elseif pool.NumWorkers ~= num_workers
    delete(pool);
    parpool('local', num_workers);
end

% Pass script directory to workers (determine before spmd)
script_dir_for_workers = script_dir;

% MPI parallel execution
spmd
    % Workers need src/ directory in their path
    % Use the script directory passed from client
    src_dir_worker = fullfile(script_dir_for_workers, 'src');
    if exist(src_dir_worker, 'dir')
        addpath(src_dir_worker);
    end
    
    % Define all constants directly inside spmd block
    % Workers cannot access client variables, so redefine everything
    halo = 1;
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
    grid = setup_simulation_grid(Np, -0.5, 0.5, -0.5, 0.5);
    M_global = setup_crossing_jets_IC(Np, Nmom_worker, rhol_worker, rhor_worker, Ma_worker, T_worker, r110_worker, r101_worker, r011_worker);
    
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
        
        % Compute fluxes for interior cells
        for i = 1:nx
            for j = 1:ny
                % Access with halo offset
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6x_hyperbolic_3D(Mr, flag2D_worker, Ma_worker);
                [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6y_hyperbolic_3D(Mr, flag2D_worker, Ma_worker);
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
        
        % Global reduction for time step (all ranks need same dt)
        vmax_local = max([abs(vpxmax(:)); abs(vpxmin(:)); abs(vpymax(:)); abs(vpymin(:))]);
        vmax = max(gcat(vmax_local, 1));  % Global max across all ranks
        dt = min(CFL_worker*dx_worker/vmax, dtmax_worker);
        dt = min(dt, tmax-t);
        t = t + dt;
        
        % X-direction flux update
        Mnpx = M;
        for j = 1:ny
            jh = j + halo;
            % Extract interior only (pas_HLL expects N×Nmom, not (N+2*halo)×Nmom)
            MOM = squeeze(M(halo+1:halo+nx, jh, :));
            FX  = squeeze(Fx(halo+1:halo+nx, jh, :));
            MNP = pas_HLL(MOM, FX, dt, dx_worker, vpxmin(:,j), vpxmax(:,j));
            Mnpx(halo+1:halo+nx, jh, :) = MNP;
        end
        
        % Y-direction flux update
        Mnpy = M;
        for i = 1:nx
            ih = i + halo;
            % Extract interior only (pas_HLL expects N×Nmom, not (N+2*halo)×Nmom)
            MOM = squeeze(M(ih, halo+1:halo+ny, :));
            FY  = squeeze(Fy(ih, halo+1:halo+ny, :));
            MNP = pas_HLL(MOM, FY, dt, dy_worker, vpymin(i,:)', vpymax(i,:)');
            Mnpy(ih, halo+1:halo+ny, :) = MNP;
        end
        
        % Combine updates (Strang splitting)
        M = Mnpx + Mnpy - M;
        
        % Exchange halos before realizability enforcement
        M = halo_exchange_2d(M, decomp, bc);
        
        % Enforce realizability
        for i = 1:nx
            for j = 1:ny
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6x_hyperbolic_3D(Mr, flag2D_worker, Ma_worker);
                [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6y_hyperbolic_3D(Mr, flag2D_worker, Ma_worker);
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
    end
    
    % Gather results to rank 1
    M_interior = M(halo+1:halo+nx, halo+1:halo+ny, :);
    
    if labindex == 1
        % Rank 1: receive from all other workers and assemble
        M_final = zeros(Np, Np, Nmom);
        M_final(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :) = M_interior;
        
        for src = 2:numlabs
            blk = labReceive(src);
            Px = decomp.dims(1);
            r = src - 1;
            rx = mod(r, Px);
            ry = floor(r / Px);
            [~, i0_s, i1_s] = block_partition_1d(Np, decomp.dims(1), rx);
            [~, j0_s, j1_s] = block_partition_1d(Np, decomp.dims(2), ry);
            M_final(i0_s:i1_s, j0_s:j1_s, :) = blk;
        end
        
        final_time = t;
        time_steps = nn;
    else
        % All other workers: send to rank 1
        labSend(M_interior, 1);
        M_final = [];
        final_time = [];
        time_steps = [];
    end
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
[C, S] = compute_CS_grid(M_final);
[M5, C5, S5] = compute_M5_grid(M_final);

% Build results structure
if nargout > 0
    results = struct();
    results.parameters.Np = Np;
    results.parameters.tmax = tmax;
    results.parameters.enable_plots = enable_plots;
    results.parameters.num_workers = num_workers;
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

function [n_local, i0, i1] = block_partition_1d(n, P, r)
    base = floor(n/P);
    remn = mod(n, P);
    if r < remn
        n_local = base + 1;
        i0 = r*(base+1) + 1;
    else
        n_local = base;
        i0 = remn*(base+1) + (r-remn)*base + 1;
    end
    i1 = i0 + n_local - 1;
end

