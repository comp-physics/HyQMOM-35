function cfg = create_config(Np, tmax, Kn, Ma, flag2D)
% create_config Create configuration struct with simulation parameters
%
% Inputs:
%   Np - number of grid points
%   tmax - maximum simulation time
%   Kn - Knudsen number
%   Ma - Mach number  
%   flag2D - 2D simulation flag
%
% Output:
%   cfg - configuration struct with all simulation parameters

cfg = struct();

% Grid parameters
cfg.Np = Np;
cfg.CFL = 0.5;
cfg.xmin = -0.5;
cfg.xmax = 0.5;
cfg.ymin = -0.5;
cfg.ymax = 0.5;

% Time parameters
cfg.tmax = tmax;
cfg.nnmax = 20000000;  % maximum number of time steps
cfg.dtmax = Kn;        % largest dt to resolve collisions

% Physical parameters
cfg.Kn = Kn;
cfg.Ma = Ma;
cfg.flag2D = flag2D;

% Moment parameters
cfg.N = 4;      % order
cfg.Nmom = 35;  % number of moments
cfg.Nmom5 = 21; % number of 5th order moments

% Initial condition parameters
cfg.rhol = 1;     % left density
cfg.rhor = 0.01;  % right density
cfg.U0 = 0;       % initial mean x-velocity
cfg.V0 = 0;       % initial mean y-velocity
cfg.W0 = 0;       % initial mean z-velocity
cfg.T = 1;        % dimensionless temperature

% Correlation coefficients for joint Gaussian
cfg.r110 = 0.0;
cfg.r101 = 0.0;
cfg.r011 = 0.0;

% Realizability parameters
cfg.s3max = 4.0 + abs(Ma)/2.0;  % limit on S300, S030, S003
cfg.h2min = 1e-8;               % minimum value for H200, H020, H002
cfg.itrealmax = 6;              % maximum realizability iterations

% Crossing parameters
cfg.Csize = floor(0.1*Np);  % size of center region

end
