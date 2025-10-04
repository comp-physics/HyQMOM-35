function varargout = grid_utils(operation, varargin)
%GRID_UTILS Unified grid and initial condition utilities
%
% Syntax:
%   grid = grid_utils('setup', Np, xmin, xmax, ymin, ymax)
%   M = grid_utils('crossing_jets_ic', Np, Nmom, rhol, rhor, Ma, T, r110, r101, r011)
%   arrays = grid_utils('init_arrays', Np, Nmom, Nmom5)
%
% Operations:
%   'setup'           - Initialize 2D spatial grid
%   'crossing_jets_ic' - Create crossing jets initial conditions
%   'init_arrays'     - Pre-allocate moment arrays

    switch operation
        case 'setup'
            varargout{1} = setup_grid_impl(varargin{:});
        case 'crossing_jets_ic'
            varargout{1} = crossing_jets_ic_impl(varargin{:});
        case 'init_arrays'
            varargout{1} = init_arrays_impl(varargin{:});
        otherwise
            error('grid_utils:unknownOp', 'Unknown operation: %s', operation);
    end
end

function grid = setup_grid_impl(Np, xmin, xmax, ymin, ymax)
% Initialize a 2D spatial grid for simulation

    grid = struct();
    
    % Generate cell edges
    grid.x = linspace(xmin, xmax, Np+1)';
    grid.y = linspace(ymin, ymax, Np+1)';
    
    % Compute cell sizes
    grid.dx = (xmax - xmin) / Np;
    grid.dy = (ymax - ymin) / Np;
    
    % Compute cell centers
    grid.xm = grid.x(1:Np) + grid.dx/2;
    grid.ym = grid.y(1:Np) + grid.dy/2;
end

function M = crossing_jets_ic_impl(Np, Nmom, rhol, rhor, Ma, T, r110, r101, r011)
% Create initial conditions for crossing jets problem

    % Initialize moment array
    M = zeros(Np, Np, Nmom);
    
    % Background parameters (mean velocities at rest)
    U0 = 0;
    V0 = 0;
    W0 = 0;
    
    % Covariance matrix for background
    C200 = T;
    C020 = T;
    C002 = T;
    C110 = r110*sqrt(C200*C020);
    C101 = r101*sqrt(C200*C002);
    C011 = r011*sqrt(C020*C002);
    
    % Initialize background (low density)
    Mr = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002);
    for i = 1:Np
        for j = 1:Np
            M(i,j,:) = Mr;
        end
    end
    
    % Crossing jets parameters
    Uc = Ma / sqrt(2);  % Velocity magnitude for crossing jets
    Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002);
    Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002);
    
    % Define jet regions (center 10% of domain on each side)
    Csize = floor(0.1*Np);
    Mint = Np/2 + 1;
    Maxt = Np/2 + 1 + Csize;
    Minb = Np/2 - Csize;
    Maxb = Np/2;
    
    % Set bottom jet (moving up-right)
    for i = Minb:Maxb
        for j = Minb:Maxb
            M(i,j,:) = Mb;
        end
    end
    
    % Set top jet (moving down-left)
    for i = Mint:Maxt
        for j = Mint:Maxt
            M(i,j,:) = Mt;
        end
    end
end

function arrays = init_arrays_impl(Np, Nmom, Nmom5)
% Pre-allocate moment arrays for simulation

    arrays = struct();
    
    % 4th-order moment arrays
    arrays.M = zeros(Np, Np, Nmom);
    arrays.C = zeros(Np, Np, Nmom);
    arrays.S = zeros(Np, Np, Nmom);
    arrays.Fx = zeros(Np, Np, Nmom);
    arrays.Fy = zeros(Np, Np, Nmom);
    
    % 5th-order moment arrays
    arrays.M5 = zeros(Np, Np, Nmom5);
    arrays.C5 = zeros(Np, Np, Nmom5);
    arrays.S5 = zeros(Np, Np, Nmom5);
    
    % Wave speed arrays
    arrays.vpxmin = zeros(Np, Np);
    arrays.vpxmax = zeros(Np, Np);
    arrays.vpymin = zeros(Np, Np);
    arrays.vpymax = zeros(Np, Np);
    arrays.v5xmin = zeros(Np, Np);
    arrays.v5xmax = zeros(Np, Np);
    arrays.v5ymin = zeros(Np, Np);
    arrays.v5ymax = zeros(Np, Np);
    arrays.v6xmin = zeros(Np, Np);
    arrays.v6xmax = zeros(Np, Np);
    arrays.v6ymin = zeros(Np, Np);
    arrays.v6ymax = zeros(Np, Np);
    
    % Temporary arrays for time stepping
    arrays.Mnp = zeros(Np, Np, Nmom);
    arrays.Mnpx = zeros(Np, Np, Nmom);
    arrays.Mnpy = zeros(Np, Np, Nmom);
    
    % Working arrays (single dimension)
    arrays.Mx = zeros(1, Nmom);
    arrays.My = zeros(1, Nmom);
    arrays.Mr = zeros(1, Nmom);
end

