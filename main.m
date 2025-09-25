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

% Handle input arguments for parameter overrides
if nargin == 0
    % Default parameters (original script behavior)
    enable_plots = true;
    save_output = false;  % Don't save by default
    Np = 6;
    tmax = 0.05;
elseif nargin == 2
    % Override Np and tmax, keep plotting enabled, no saving
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = true;
    save_output = false;
elseif nargin == 3
    % Override all three main parameters, no saving
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = varargin{3};
    save_output = false;
elseif nargin == 4
    % Override all parameters including save option
    Np = varargin{1};
    tmax = varargin{2};
    enable_plots = varargin{3};
    save_output = varargin{4};
else
    error('Invalid number of arguments. Usage: main() or main(Np, tmax) or main(Np, tmax, enable_plots) or main(Np, tmax, enable_plots, save_output)');
end

% Fixed simulation parameters
% Knudsen number (>= 0.001 to avoid long simulations)
Kn = 1/1;

% Mach number (for impinging jets with velocity u and temperature Theta)
Ma = 0;  % (= u/sqrt(Theta))

% flag for 2-D case (use only if S101=S011=0) if flag2D == 1
flag2D = 0;

%% 2-D space discretization: square domain
CFL = 0.5;
xmin = -0.5;
xmax = 0.5;
ymin = -0.5;
ymax = 0.5;
x = xmin + (xmax-xmin)*linspace(0,1,Np+1)';
y = ymin + (ymax-ymin)*linspace(0,1,Np+1)';
dx = (xmax-xmin)/Np;
xm = x(1:Np)+dx/2;
dy = (ymax-ymin)/Np;
ym = y(1:Np)+dy/2;
%%

% order and number of moments (fixed)
N = 4;
Nmom = 35;
Nmom5 = 21;

% moment index map to avoid magic numbers
idx = moment_indices();

% maximum number of time steps
nnmax = 20000000;
%nnmax = 5;

% initial correlation coefficients for joint Gaussian
r110 = 0.;
r101 = 0.;
r011 = 0.;

% largest dt to resolve collisions
dtmax = Kn;

%% shock problem %%%%%%%%%%%%
% initial densities
rhol = 1;
rhor = 0.01;

% initial mean velocities
U0 = 0;
V0 = 0;
W0 = 0;

% dimensionless temperature: used for scaling velocity so T=1 (do not change)
T = 1;

% set initial conditions to joint Gaussian with covariance
C200 = T;
C020 = T;
C002 = T;
C110 = r110*sqrt(C200*C020);
C101 = r101*sqrt(C200*C002);
C011 = r011*sqrt(C020*C002);

% initialize moments on "left" and "right"
Ml = InitializeM4_35(rhol,U0,V0,W0,C200,C110,C101,C020,C011,C002);
Mr = InitializeM4_35(rhor,U0,V0,W0,C200,C110,C101,C020,C011,C002);

%%
C200c = T;
C020c = T;
C002c = T;
C110c = r110*sqrt(C200*C020);
C101c = r101*sqrt(C200*C002);
C011c = r011*sqrt(C020*C002);
% magnitude of 3-D velocity = Ma
Uc = Ma/sqrt(2);
% initialize moments on "top" and "bottom" for crossing
Mt = InitializeM4_35(rhol,-Uc,-Uc,W0,C200c,C110c,C101c,C020c,C011c,C002c);
Mb = InitializeM4_35(rhol, Uc, Uc,W0,C200c,C110c,C101c,C020c,C011c,C002c);
%%

M = zeros(Np,Np,Nmom);
Mnp = M;
Mnpx = M;
Mnpy = M;
S = zeros(Np,Np,Nmom);
C = zeros(Np,Np,Nmom);
M5 = zeros(Np,Np,Nmom5);
S5 = zeros(Np,Np,Nmom5);
C5 = zeros(Np,Np,Nmom5);

%% initialize 35 3-D moments on 2-D spatial domain
% low-pressure background
M = repmat(reshape(Mr,1,1,[]), Np, Np, 1);

% high-pressure center (Csize = size of center region)
Csize = floor(0.1*Np) ;
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
txt = ['riemann_3D_hyqmom35_crossing','_Np',num2str(Np),'_Kn',num2str(Kn),'_Ma',num2str(Ma),'.mat'];

%% time evolution begins here
t = 0.;
Fx = zeros(Np,Np,Nmom);
Fy = zeros(Np,Np,Nmom);
Mx = zeros(1,Nmom);  % closures for x flux
My = zeros(1,Nmom);  % closures for y flux
Mr = zeros(1,Nmom);  % realizable moments
vpxmin = zeros(Np,Np,1);
vpxmax = zeros(Np,Np,1);
vpymin = zeros(Np,Np,1);
vpymax = zeros(Np,Np,1);
v5xmin = zeros(Np,Np,1);
v5xmax = zeros(Np,Np,1);
v5ymin = zeros(Np,Np,1);
v5ymax = zeros(Np,Np,1);
v6xmin = zeros(Np,Np,1);
v6xmax = zeros(Np,Np,1);
v6ymin = zeros(Np,Np,1);
v6ymax = zeros(Np,Np,1);
nn = 0;

tic
while t<tmax && nn<nnmax
    nn = nn+1;
    
    % spatial fluxes, realizability checks, eigenvalues
    Mnp = M;
    parfor i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i,j,:));
            % eigenvalues with hyperbolicity
            [Mx,My,~,Mr] = Flux_closure35_and_realizable_3D(MOM,flag2D,Ma);
            [v6xmin(i,j),v6xmax(i,j),Mr] = eigenvalues6_hyperbolic_3D(Mr,'x',flag2D,Ma);
            [v6ymin(i,j),v6ymax(i,j),Mr] = eigenvalues6_hyperbolic_3D(Mr,'y',flag2D,Ma);
            % fluxes in the x direction
            Fx(i,j,:) = Mx;
            % fluxes in the y direction
            Fy(i,j,:) = My;
            % realizable moments
            Mnp(i,j,:)= Mr;
            %
            % compute eigenvalues for HLL
            % 1-D hyqmom for m500 eigenvalues in x direction
            MOM5 = Mr(idx.x_moments); % m000 m100 m200 m300 m400
            [~,v5xmin(i,j),v5xmax(i,j)] = closure_and_eigenvalues(MOM5);
            %
            vpxmin(i,j)=min(v5xmin(i,j),v6xmin(i,j));
            vpxmax(i,j)=max(v5xmax(i,j),v6xmax(i,j));
            % 1-D hyqmom for m050 eigenvalues in y direction
            MOM5 = Mr(idx.y_moments); % m000 m010 m020 m030 m040
            [~,v5ymin(i,j),v5ymax(i,j)] = closure_and_eigenvalues(MOM5);
            %
            vpymin(i,j)=min(v5ymin(i,j),v6ymin(i,j));
            vpymax(i,j)=max(v5ymax(i,j),v6ymax(i,j));
            % 
        end
    end
    M = Mnp;

    % fix time step based on largest eigenvalues in computational domain
    dt = CFL*dx/max([abs(vpxmax);abs(vpxmin);abs(vpymax);abs(vpymin)],[],'all');
    if t+dt>tmax
        dt=tmax-t;
    end
    dt = min(dt,dtmax);
    t = t+dt;
    
    %% Euler for flux starts here
    % update moments due to spatial fluxes using method of lines and HLL
    Mnpx = apply_hll_update(M, Fx, vpxmin, vpxmax, dt, dx, 'x');
    Mnpy = apply_hll_update(M, Fy, vpymin, vpymax, dt, dy, 'y');
    % end of Euler (NB: Mnp can have unrealizable moments)
    Mnp = Mnpx + Mnpy - M;
    %%
    M = Mnp;
    %
    % enforce realizability and hyperbolicity
    parfor i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i,j,:));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM,flag2D,Ma);
            [v6xmin(i,j),v6xmax(i,j),Mr] = eigenvalues6_hyperbolic_3D(Mr,'x',flag2D,Ma);
            [v6ymin(i,j),v6ymax(i,j),Mr] = eigenvalues6_hyperbolic_3D(Mr,'y',flag2D,Ma);
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr,flag2D,Ma);
            % realizable moments
            Mnp(i,j,:)= Mr;
        end
    end
    M = Mnp;
    %
    % collision step using BGK
    parfor i = 1:Np
        for j = 1:Np
            MM = squeeze(M(i,j,:));
            MMC = collision35(MM,dt,Kn);
            Mnp(i,j,:) = MMC;
        end
    end
    M = Mnp;
    %
    % compute central and standardized moments, check 1-D realizability
    [C, S] = compute_CS_grid(M);
    %
    if any(C(:,:,idx.C200) < 0,'all') 
        disp('pb C200 realizabilite apres pas temps')
        break
    end
    if any(C(:,:,idx.C020) < 0,'all')
        disp('pb C020 realizabilite apres pas temps')
        break
    end
    if any(C(:,:,idx.C002) < 0,'all')
        disp('pb C002 realizabilite apres pas temps')
        break
    end
    %
    if any(S(:,:,5)-1-S(:,:,4).^2 < 0,'all') 
        disp('pb H200 realizabilite apres pas temps')
        %break
    end
    if any(S(:,:,15)-1-S(:,:,13).^2 < 0,'all') 
        disp('pb H020 realizabilite apres pas temps')
        %break
    end
    if any(S(:,:,25)-1-S(:,:,23).^2 < 0,'all') 
        disp('pb H002 realizabilite apres pas temps')
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