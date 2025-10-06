function params = parse_main_args(args, defaults)
% Parse command-line arguments for main simulation
% Uses inputParser for clean, extensible argument handling
%
% Inputs:
%   args     - Cell array of input arguments (varargin from main)
%   defaults - Struct with default values for all parameters
%
% Output:
%   params   - Struct with all parsed parameters

p = inputParser;
p.FunctionName = 'main';

% Add all parameters with defaults
addOptional(p, 'Np', defaults.Np, @(x) isnumeric(x) && isscalar(x) && x > 0);
addOptional(p, 'tmax', defaults.tmax, @(x) isnumeric(x) && isscalar(x) && x > 0);
addOptional(p, 'enable_plots', defaults.enable_plots, @islogical);
addOptional(p, 'num_workers', defaults.num_workers, @(x) isnumeric(x) && isscalar(x) && x > 0);
addOptional(p, 'enable_profile', defaults.enable_profile, @islogical);
addOptional(p, 'symmetry_check_interval', defaults.symmetry_check_interval, @(x) isnumeric(x) && isscalar(x) && x > 0);

% Parse input
parse(p, args{:});

% Extract to struct
params = p.Results;

% Add physical parameters (constants)
params.Kn = 1.0;
params.Ma = 0.0;
params.flag2D = 0;
params.CFL = 0.5;
params.dx = 1.0 / params.Np;
params.dy = 1.0 / params.Np;
params.N = 4;
params.Nmom = 35;
params.Nmom5 = 21;
params.nnmax = 2e7;
params.dtmax = params.Kn;

% Correlation coefficients
params.r110 = 0.0;
params.r101 = 0.0;
params.r011 = 0.0;

% Initial conditions
params.T = 1.0;
params.rhol = 1.0;
params.rhor = 0.01;

end
