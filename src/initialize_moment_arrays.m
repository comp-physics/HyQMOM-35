function arrays = initialize_moment_arrays(Np, Nmom, Nmom5)
% INITIALIZE_MOMENT_ARRAYS Pre-allocates moment arrays for simulation
%
%   arrays = initialize_moment_arrays(Np, Nmom, Nmom5)
%
%   Inputs:
%       Np - Grid points per dimension
%       Nmom - Number of 4th-order moments (35)
%       Nmom5 - Number of 5th-order moments (21)
%
%   Output:
%       arrays - Structure with pre-allocated arrays

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

