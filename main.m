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

% Parse input arguments with defaults (before any clearing)
p = inputParser;
addOptional(p, 'Np', 6, @(x) isnumeric(x) && isscalar(x) && x > 0);
addOptional(p, 'tmax', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
addOptional(p, 'enable_plots', false, @(x) islogical(x) || isnumeric(x));
addOptional(p, 'save_output', false, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

% Extract parsed parameters
Np = p.Results.Np;
tmax = p.Results.tmax;
enable_plots = logical(p.Results.enable_plots);
save_output = logical(p.Results.save_output);

% Clear and initialize (only if running as script)
if nargin == 0
    clc
    % Clear all variables except the parsed parameters
    clearvars('-except', 'Np', 'tmax', 'enable_plots', 'save_output')
    close all
end

% Add src directory to path for function dependencies
% Get the directory where this script is located
script_dir = fileparts(mfilename('fullpath'));
src_dir = fullfile(script_dir, 'src');
autogen_dir = fullfile(src_dir, 'autogen');
if exist(src_dir, 'dir')
    addpath(src_dir);
end
if exist(autogen_dir, 'dir')
    addpath(autogen_dir);
end

% Create configuration struct with all simulation parameters
cfg = create_config(Np, tmax, 1.0, 0.0, 0);  % Kn=1, Ma=0, flag2D=0

% Configuration struct contains all parameters - use cfg.* directly

%% 2-D space discretization: square domain
x = cfg.xmin + (cfg.xmax-cfg.xmin)*linspace(0,1,Np+1)';
y = cfg.ymin + (cfg.ymax-cfg.ymin)*linspace(0,1,Np+1)';
dx = (cfg.xmax-cfg.xmin)/Np;
xm = x(1:Np)+dx/2;
dy = (cfg.ymax-cfg.ymin)/Np;
ym = y(1:Np)+dy/2;

% Use cfg.N, cfg.Nmom, cfg.Nmom5 directly

% moment index map to avoid magic numbers
idx = moment_indices();

%% shock problem %%%%%%%%%%%%
% Use cfg initial condition parameters directly

% set initial conditions to joint Gaussian with covariance
C200 = cfg.T;
C020 = cfg.T;
C002 = cfg.T;
C110 = cfg.r110*sqrt(C200*C020);
C101 = cfg.r101*sqrt(C200*C002);
C011 = cfg.r011*sqrt(C020*C002);

% initialize moments on "left" and "right"
Ml = InitializeM4_35(cfg.rhol,cfg.U0,cfg.V0,cfg.W0,C200,C110,C101,C020,C011,C002);
Mr = InitializeM4_35(cfg.rhor,cfg.U0,cfg.V0,cfg.W0,C200,C110,C101,C020,C011,C002);

%%
C200c = cfg.T;
C020c = cfg.T;
C002c = cfg.T;
C110c = cfg.r110*sqrt(C200*C020);
C101c = cfg.r101*sqrt(C200*C002);
C011c = cfg.r011*sqrt(C020*C002);
% magnitude of 3-D velocity = Ma
Uc = cfg.Ma/sqrt(2);
% initialize moments on "top" and "bottom" for crossing
Mt = InitializeM4_35(cfg.rhol,-Uc,-Uc,cfg.W0,C200c,C110c,C101c,C020c,C011c,C002c);
Mb = InitializeM4_35(cfg.rhol, Uc, Uc,cfg.W0,C200c,C110c,C101c,C020c,C011c,C002c);
%%

M = zeros(Np,Np,cfg.Nmom);
Mnp = M;
Mnpx = M;
Mnpy = M;
S = zeros(Np,Np,cfg.Nmom);
C = zeros(Np,Np,cfg.Nmom);
M5 = zeros(Np,Np,cfg.Nmom5);
S5 = zeros(Np,Np,cfg.Nmom5);
C5 = zeros(Np,Np,cfg.Nmom5);

%% initialize 35 3-D moments on 2-D spatial domain
% low-pressure background
M = repmat(reshape(Mr,1,1,[]), Np, Np, 1);

% high-pressure center (Csize = size of center region)
Csize = cfg.Csize;
Mint = Np/2 + 1;
Maxt = Np/2 + 1 + Csize;
Minb = Np/2 - Csize;
Maxb = Np/2;
M(Minb:Maxb, Minb:Maxb, :) = repmat(reshape(Mb,1,1,[]), Maxb-Minb+1, Maxb-Minb+1, 1);
M(Mint:Maxt, Mint:Maxt, :) = repmat(reshape(Mt,1,1,[]), Maxt-Mint+1, Maxt-Mint+1, 1);
%
% compute moments and plot initial conditions
[C, S] = compute_CS_grid(M);
[M5, C5, S5] = compute_M5_grid(M);

nmin = 1;
nmax = Np;

cc = 'k';

% Plot initial conditions
simulation_plots('initial', xm, ym, M, C, S, M5, C5, S5, Np, enable_plots);
%%

% name saved file
txt = ['riemann_3D_hyqmom35_crossing','_Np',num2str(Np),'_Kn',num2str(cfg.Kn),'_Ma',num2str(cfg.Ma),'.mat'];

%% time evolution begins here
t = 0.;
Fx = zeros(Np,Np,cfg.Nmom);
Fy = zeros(Np,Np,cfg.Nmom);
Mr = zeros(1,cfg.Nmom);  % realizable moments

% Consolidated bounds storage using struct arrays
bounds_grid = struct('hll', struct('xmin', zeros(Np,Np), 'xmax', zeros(Np,Np), ...
                                   'ymin', zeros(Np,Np), 'ymax', zeros(Np,Np)), ...
                     'x', struct('v6min', zeros(Np,Np), 'v6max', zeros(Np,Np), ...
                                 'v5min', zeros(Np,Np), 'v5max', zeros(Np,Np)), ...
                     'y', struct('v6min', zeros(Np,Np), 'v6max', zeros(Np,Np), ...
                                 'v5min', zeros(Np,Np), 'v5max', zeros(Np,Np)));
nn = 0;

tic
while t<tmax && nn<cfg.nnmax
    nn = nn+1;
    
    % spatial fluxes, realizability checks, eigenvalues
    Mnp = M;
    
    % Extract bounds arrays for parfor compatibility
    v6xmin = bounds_grid.x.v6min;
    v6xmax = bounds_grid.x.v6max;
    v5xmin = bounds_grid.x.v5min;
    v5xmax = bounds_grid.x.v5max;
    v6ymin = bounds_grid.y.v6min;
    v6ymax = bounds_grid.y.v6max;
    v5ymin = bounds_grid.y.v5min;
    v5ymax = bounds_grid.y.v5max;
    vpxmin = bounds_grid.hll.xmin;
    vpxmax = bounds_grid.hll.xmax;
    vpymin = bounds_grid.hll.ymin;
    vpymax = bounds_grid.hll.ymax;
    
    parfor i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i,j,:));
            [Mr, flux, bounds] = process_cell_timestep(MOM, cfg.flag2D, cfg.Ma, idx);
            
            % Store results using structured format
            Fx(i,j,:) = flux.x;
            Fy(i,j,:) = flux.y;
            Mnp(i,j,:) = Mr;
            
            % Store bounds in temporary arrays
            v6xmin(i,j) = bounds.x.v6min;
            v6xmax(i,j) = bounds.x.v6max;
            v5xmin(i,j) = bounds.x.v5min;
            v5xmax(i,j) = bounds.x.v5max;
            v6ymin(i,j) = bounds.y.v6min;
            v6ymax(i,j) = bounds.y.v6max;
            v5ymin(i,j) = bounds.y.v5min;
            v5ymax(i,j) = bounds.y.v5max;
            vpxmin(i,j) = bounds.hll.xmin;
            vpxmax(i,j) = bounds.hll.xmax;
            vpymin(i,j) = bounds.hll.ymin;
            vpymax(i,j) = bounds.hll.ymax;
        end
    end
    
    % Update bounds_grid structure after parfor
    bounds_grid.x.v6min = v6xmin;
    bounds_grid.x.v6max = v6xmax;
    bounds_grid.x.v5min = v5xmin;
    bounds_grid.x.v5max = v5xmax;
    bounds_grid.y.v6min = v6ymin;
    bounds_grid.y.v6max = v6ymax;
    bounds_grid.y.v5min = v5ymin;
    bounds_grid.y.v5max = v5ymax;
    bounds_grid.hll.xmin = vpxmin;
    bounds_grid.hll.xmax = vpxmax;
    bounds_grid.hll.ymin = vpymin;
    bounds_grid.hll.ymax = vpymax;
    
    M = Mnp;

    % fix time step based on largest eigenvalues in computational domain
    dt = cfg.CFL*dx/max([abs(bounds_grid.hll.xmax(:)); abs(bounds_grid.hll.xmin(:)); ...
                     abs(bounds_grid.hll.ymax(:)); abs(bounds_grid.hll.ymin(:))]);
    if t+dt>tmax
        dt=tmax-t;
    end
    dt = min(dt,cfg.dtmax);
    t = t+dt;
    
    %% Euler for flux starts here
    % update moments due to spatial fluxes using method of lines and HLL
    Mnpx = apply_hll_update(M, Fx, bounds_grid.hll.xmin, bounds_grid.hll.xmax, dt, dx, 'x');
    Mnpy = apply_hll_update(M, Fy, bounds_grid.hll.ymin, bounds_grid.hll.ymax, dt, dy, 'y');
    % end of Euler (NB: Mnp can have unrealizable moments)
    Mnp = Mnpx + Mnpy - M;
    %%
    M = Mnp;
    %
    % enforce realizability and hyperbolicity
    v6xmin = bounds_grid.x.v6min;
    v6xmax = bounds_grid.x.v6max;
    v6ymin = bounds_grid.y.v6min;
    v6ymax = bounds_grid.y.v6max;
    
    parfor i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i,j,:));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM,cfg.flag2D,cfg.Ma);
            [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr,'x',cfg.flag2D,cfg.Ma);
            [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr,'y',cfg.flag2D,cfg.Ma);
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr,cfg.flag2D,cfg.Ma);
            % realizable moments
            Mnp(i,j,:)= Mr;
        end
    end
    
    % Update bounds_grid structure after parfor
    bounds_grid.x.v6min = v6xmin;
    bounds_grid.x.v6max = v6xmax;
    bounds_grid.y.v6min = v6ymin;
    bounds_grid.y.v6max = v6ymax;
    
    M = Mnp;
    %
    % collision step using BGK
    parfor i = 1:Np
        for j = 1:Np
            MM = squeeze(M(i,j,:));
            MMC = collision35(MM,dt,cfg.Kn);
            Mnp(i,j,:) = MMC;
        end
    end
    M = Mnp;
    %
    % compute central and standardized moments, check 1-D realizability
    [C, S] = compute_CS_grid(M);
    %
    if any(C(:,:,idx.C200) < 0,'all') 
        warning('C200 < 0 after timestep %d at t=%.6f; aborting simulation', nn, t);
        break
    end
    if any(C(:,:,idx.C020) < 0,'all')
        warning('C020 < 0 after timestep %d at t=%.6f; aborting simulation', nn, t);
        break
    end
    if any(C(:,:,idx.C002) < 0,'all')
        warning('C002 < 0 after timestep %d at t=%.6f; aborting simulation', nn, t);
        break
    end
    %
    if any(S(:,:,5)-1-S(:,:,4).^2 < 0,'all') 
        warning('H200 realizability violation after timestep %d at t=%.6f', nn, t);
        %break
    end
    if any(S(:,:,15)-1-S(:,:,13).^2 < 0,'all') 
        warning('H020 realizability violation after timestep %d at t=%.6f', nn, t);
        %break
    end
    if any(S(:,:,25)-1-S(:,:,23).^2 < 0,'all') 
        warning('H002 realizability violation after timestep %d at t=%.6f', nn, t);
        %break
    end

    % Plot time evolution
    simulation_plots('time_evolution', S, C, xm, ym, Np, enable_plots);

    [Diff,MaxDiff] = test_symmetry_2D(M,Np);
    MaxDiff 
end
toc

%% postprocessing for plots
[M5, C5, S5] = compute_M5_grid(M);

[lam6xa, lam6xb, lam6ya, lam6yb] = compute_jacobian_eigenvalues(M);

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
simulation_plots('final', xm, ym, M, C, S, M5, C5, S5, Np, ...
                 bounds_grid.x.v5min, bounds_grid.x.v5max, bounds_grid.x.v6min, bounds_grid.x.v6max, ...
                 bounds_grid.y.v5min, bounds_grid.y.v5max, bounds_grid.y.v6min, bounds_grid.y.v6max, ...
                 lam6xa, lam6xb, lam6ya, lam6yb, enable_plots);

% Return results structure (only if output is requested)
if nargout > 0
    results = struct();
    
    % Simulation parameters
    results.parameters.Np = Np;
    results.parameters.tmax = tmax;
    results.parameters.enable_plots = enable_plots;
    results.parameters.save_output = save_output;
    results.parameters.Kn = cfg.Kn;
    results.parameters.Ma = cfg.Ma;
    results.parameters.CFL = cfg.CFL;
    results.parameters.Nmom = cfg.Nmom;
    results.parameters.N = cfg.N;
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
    
    % Velocity bounds (consolidated structure)
    results.bounds = bounds_grid;
    
    % Filename for saving
    results.filename = txt;
end

end