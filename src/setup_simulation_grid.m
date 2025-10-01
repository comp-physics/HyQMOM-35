function grid = setup_simulation_grid(Np, xmin, xmax, ymin, ymax)
% SETUP_SIMULATION_GRID Initializes a 2D spatial grid for simulation
%
%   grid = setup_simulation_grid(Np, xmin, xmax, ymin, ymax)
%
%   Inputs:
%       Np - Number of grid points per dimension
%       xmin, xmax - x-direction domain bounds
%       ymin, ymax - y-direction domain bounds
%
%   Output:
%       grid - Structure containing grid information:
%           .x, .y   - Cell edge coordinates
%           .xm, .ym - Cell center coordinates
%           .dx, .dy - Cell sizes

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

