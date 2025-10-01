function [results] = main(varargin)
% 2-D Riemann solver for 3-D HyQMOM using HLL and explicit Euler
% code restricted to N=4 in 3-D with 35 moments
%
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]
%
% Usage:
%   main()                    % Run with default parameters
%   main(Np, tmax)           % Override Np and tmax
%   main(Np, tmax, enable_plots) % Override plotting
%   main(Np, tmax, enable_plots, save_output) % Override all parameters
%
% Examples:
%   main()           % Default: Np=6, tmax=0.05, enable_plots=true, save_output=false
%   main(6, 0.02)    % Golden file parameters with plotting enabled, no saving
%   main(6, 0.02, false) % Golden file parameters without plotting, no saving
%   main(6, 0.02, false, true) % Golden file creation with saving enabled
%
% This version limits the energy fluxes and checks realizability and
% hyperbolicity

% Clear and initialize (only if running as script)
if nargin == 0
    clc
    clear 
    close all
end

% Add src directory to path for function dependencies
% Get the directory where this script is located
script_dir = fileparts(mfilename('fullpath'));
src_dir = fullfile(script_dir, 'src');
if exist(src_dir, 'dir')
    addpath(src_dir);
end

% Parse input arguments with defaults
defaults = struct('Np', 10, 'tmax', 0.1, 'enable_plots', false, 'save_output', false);
if nargin == 0
    Np = defaults.Np;
    tmax = defaults.tmax;
    enable_plots = defaults.enable_plots;
    save_output = defaults.save_output;
else
    [Np, tmax, enable_plots, save_output] = parse_input_args(nargin, varargin, defaults);
end

% Physical parameters
Kn = 1.0;      % Knudsen number
Ma = 0.0;      % Mach number
flag2D = 0;    % 2D mode flag (0 = full 3D)

%% Spatial and temporal discretization
CFL = 0.5;
grid = setup_simulation_grid(Np, -0.5, 0.5, -0.5, 0.5);
x = grid.x; y = grid.y; xm = grid.xm; ym = grid.ym; dx = grid.dx; dy = grid.dy;

N = 4;          % Moment order
Nmom = 35;      % Number of 4th-order moments
Nmom5 = 21;     % Number of 5th-order moments
nnmax = 2e7;    % Maximum time steps
dtmax = Kn;     % Maximum time step for collision resolution

% Correlation coefficients for initial Gaussian
r110 = 0.0;
r101 = 0.0;
r011 = 0.0;

%% Initial conditions: Crossing jets problem
T = 1.0;        % Dimensionless temperature
rhol = 1.0;     % High density (jets)
rhor = 0.01;    % Low density (background)

% Initialize moment arrays
arrays = initialize_moment_arrays(Np, Nmom, Nmom5);
M = setup_crossing_jets_IC(Np, Nmom, rhol, rhor, Ma, T, r110, r101, r011);

% Compute derived moments
[C, S] = compute_CS_grid(M);
[M5, C5, S5] = compute_M5_grid(M);

% Extract pre-allocated arrays
Mnp = arrays.Mnp; Mnpx = arrays.Mnpx; Mnpy = arrays.Mnpy;
Fx = arrays.Fx; Fy = arrays.Fy;
vpxmin = arrays.vpxmin; vpxmax = arrays.vpxmax;
vpymin = arrays.vpymin; vpymax = arrays.vpymax;
v5xmin = arrays.v5xmin; v5xmax = arrays.v5xmax;
v5ymin = arrays.v5ymin; v5ymax = arrays.v5ymax;
v6xmin = arrays.v6xmin; v6xmax = arrays.v6xmax;
v6ymin = arrays.v6ymin; v6ymax = arrays.v6ymax;
Mx = arrays.Mx; My = arrays.My; Mr = arrays.Mr;

nmin = 1;
nmax = Np;

cc = 'k';

% Plot initial conditions
simulation_plots('initial', xm, ym, M, C, S, M5, C5, S5, Np, enable_plots);
%%

% Output filename
txt = sprintf('riemann_3D_hyqmom35_crossing_Np%d_Kn%g_Ma%g.mat', Np, Kn, Ma);

%% Time evolution
t = 0.0;
nn = 0;

tic
while t<tmax && nn<nnmax
    nn = nn+1;
    
    % Compute fluxes, realizability, and eigenvalues for each grid point
    Mnp = M;
    for i = 1:Np  % Changed from parfor: MPI will handle parallelism
        for j = 1:Np
            MOM = squeeze(M(i,j,:));
            
            % Enforce realizability and compute fluxes
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6x_hyperbolic_3D(Mr, flag2D, Ma);
            [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6y_hyperbolic_3D(Mr, flag2D, Ma);
            [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            
            Fx(i,j,:) = Mx;
            Fy(i,j,:) = My;
            Mnp(i,j,:) = Mr;
            
            % Compute 1D eigenvalue bounds for HLL solver
            [~, v5xmin(i,j), v5xmax(i,j)] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
            [~, v5ymin(i,j), v5ymax(i,j)] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
            
            % Combined eigenvalue bounds
            vpxmin(i,j) = min(v5xmin(i,j), v6xmin(i,j));
            vpxmax(i,j) = max(v5xmax(i,j), v6xmax(i,j));
            vpymin(i,j) = min(v5ymin(i,j), v6ymin(i,j));
            vpymax(i,j) = max(v5ymax(i,j), v6ymax(i,j));
        end
    end
    M = Mnp;

    % Compute time step from CFL condition
    vmax = max([abs(vpxmax(:)); abs(vpxmin(:)); abs(vpymax(:)); abs(vpymin(:))]);
    dt = min(CFL*dx/vmax, dtmax);
    dt = min(dt, tmax-t);  % Don't overshoot tmax
    t = t + dt
    
    % X-direction flux update using dimensional splitting
    Mnpx = M;
    for j = 1:Np  % Changed from parfor: MPI will handle parallelism
        MOM = squeeze(M(:,j,:));
        FX = squeeze(Fx(:,j,:));
        MNP = pas_HLL(MOM, FX, dt, dx, vpxmin(:,j), vpxmax(:,j));
        Mnpx(:,j,:) = MNP;
    end
    
    % Y-direction flux update
    Mnpy = M;
    for i = 1:Np  % Changed from parfor: MPI will handle parallelism
        MOM = squeeze(M(i,:,:));
        FY = squeeze(Fy(i,:,:));
        MNP = pas_HLL(MOM, FY, dt, dy, vpymin(i,:)', vpymax(i,:)');
        Mnpy(i,:,:) = MNP;
    end
    
    % Combine updates (Strang splitting: Mnp = Lx(Ly(M)))
    M = Mnpx + Mnpy - M;
    % Enforce realizability and hyperbolicity after advection
    for i = 1:Np  % Changed from parfor: MPI will handle parallelism
        for j = 1:Np
            MOM = squeeze(M(i,j,:));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6x_hyperbolic_3D(Mr, flag2D, Ma);
            [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6y_hyperbolic_3D(Mr, flag2D, Ma);
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            Mnp(i,j,:) = Mr;
        end
    end
    M = Mnp;
    
    % Apply BGK collision operator
    for i = 1:Np  % Changed from parfor: MPI will handle parallelism
        for j = 1:Np
            MOM = squeeze(M(i,j,:));
            MMC = collision35(MOM, dt, Kn);
            Mnp(i,j,:) = MMC;
        end
    end
    M = Mnp;
    % Check realizability constraints
    [C, S] = compute_CS_grid(M);
    
    % Check variances (must be positive)
    if any(C(:,:,3) < 0, 'all'), error('C200 realizability violation'); end
    if any(C(:,:,10) < 0, 'all'), error('C020 realizability violation'); end
    if any(C(:,:,20) < 0, 'all'), error('C002 realizability violation'); end
    
    % Check hyperbolic realizability constraints (kurtosis bounds)
    if any(S(:,:,5) - 1 - S(:,:,4).^2 < 0, 'all')
        warning('H200 realizability constraint violated');
    end
    if any(S(:,:,15) - 1 - S(:,:,13).^2 < 0, 'all')
        warning('H020 realizability constraint violated');
    end
    if any(S(:,:,25) - 1 - S(:,:,23).^2 < 0, 'all')
        warning('H002 realizability constraint violated');
    end

    % Plot time evolution
    simulation_plots('time_evolution', S, C, xm, ym, Np, enable_plots);

    [Diff,MaxDiff] = test_symmetry_2D(M,Np);
    MaxDiff 
end
toc

%% postprocessing for plots
[M5, C5, S5] = compute_M5_grid(M);

% Compute Jacobian eigenvalues for postprocessing plots
lam6xa = zeros(Np,Np,6);
lam6xb = zeros(Np,Np,6);
lam6ya = zeros(Np,Np,6);
lam6yb = zeros(Np,Np,6);
for i = 1:Np
    for j = 1:Np
        M1 = squeeze(M(i,j,:));
        % X-direction eigenvalues (UV and UW planes)
        moments_uv = [M1(1),M1(6),M1(10),M1(13),M1(15),M1(2),M1(7),M1(11),M1(14),M1(3),M1(8),M1(12),M1(4),M1(9),M1(5)];
        moments_uw = [M1(1),M1(16),M1(20),M1(23),M1(25),M1(2),M1(17),M1(21),M1(24),M1(3),M1(18),M1(22),M1(4),M1(19),M1(5)];
        [~, ~, lam6xa_temp, lam6xb_temp] = compute_jacobian_eigenvalues(moments_uv, moments_uw);
        lam6xa(i,j,:) = lam6xa_temp;
        lam6xb(i,j,:) = lam6xb_temp;
        % Y-direction eigenvalues (VU and VW planes)
        moments_vu = [M1(1),M1(2),M1(3),M1(4),M1(5),M1(6),M1(7),M1(8),M1(9),M1(10),M1(11),M1(12),M1(13),M1(14),M1(15)];
        moments_vw = [M1(1),M1(16),M1(20),M1(23),M1(25),M1(6),M1(26),M1(32),M1(34),M1(10),M1(29),M1(35),M1(13),M1(31),M1(15)];
        [~, ~, lam6ya_temp, lam6yb_temp] = compute_jacobian_eigenvalues(moments_vu, moments_vw);
        lam6ya(i,j,:) = lam6ya_temp;
        lam6yb(i,j,:) = lam6yb_temp;
    end
end

% save simulation data (only if requested)
if save_output
    save(txt)
    fprintf('Simulation data saved to: %s\n', txt);
end
%%

%% plots
nmin = 1;
nmax = Np;

cc = 'r';

% Plot final results
simulation_plots('final', xm, ym, M, C, S, M5, C5, S5, Np, v5xmin, v5xmax, v6xmin, v6xmax, v5ymin, v5ymax, v6ymin, v6ymax, lam6xa, lam6xb, lam6ya, lam6yb, enable_plots);

% Return results structure (only if output is requested)
if nargout > 0
    results = struct();
    
    % Simulation parameters
    results.parameters.Np = Np;
    results.parameters.tmax = tmax;
    results.parameters.enable_plots = enable_plots;
    results.parameters.save_output = save_output;
    results.parameters.Kn = Kn;
    results.parameters.Ma = Ma;
    results.parameters.CFL = CFL;
    results.parameters.Nmom = Nmom;
    results.parameters.N = N;
    results.parameters.final_time = t;
    results.parameters.time_steps = nn;
    
    % Spatial grid
    results.grid.x = x;
    results.grid.y = y;
    results.grid.xm = xm;
    results.grid.ym = ym;
    results.grid.dx = dx;
    results.grid.dy = dy;
    
    % Final moment data
    results.moments.M = M;
    results.moments.C = C;
    results.moments.S = S;
    results.moments.M5 = M5;
    results.moments.C5 = C5;
    results.moments.S5 = S5;
    
    % Eigenvalue data
    if exist('lam6xa', 'var')
        results.eigenvalues.lam6xa = lam6xa;
        results.eigenvalues.lam6xb = lam6xb;
        results.eigenvalues.lam6ya = lam6ya;
        results.eigenvalues.lam6yb = lam6yb;
    end
    
    % Velocity bounds
    if exist('v5xmin', 'var')
        results.velocities.v5xmin = v5xmin;
        results.velocities.v5xmax = v5xmax;
        results.velocities.v5ymin = v5ymin;
        results.velocities.v5ymax = v5ymax;
        results.velocities.v6xmin = v6xmin;
        results.velocities.v6xmax = v6xmax;
        results.velocities.v6ymin = v6ymin;
        results.velocities.v6ymax = v6ymax;
        results.velocities.vpxmin = vpxmin;
        results.velocities.vpxmax = vpxmax;
        results.velocities.vpymin = vpymin;
        results.velocities.vpymax = vpymax;
    end
    
    % Filename for saving
    results.filename = txt;
end

end