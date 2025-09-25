function [Mr_final, flux, bounds] = process_cell_timestep(MOM, flag2D, Ma, idx)
% process_cell_timestep Consolidated per-cell operations for timestep
%
% This function combines the per-cell operations for flux computation and
% eigenvalue calculation only. Realizability enforcement and collision
% are handled separately to maintain the original order of operations.
%
% Inputs:
%   MOM - moment vector for this cell (Nmom x 1)
%   flag2D - 2D flag
%   Ma - Mach number
%   idx - moment indices structure
%
% Outputs:
%   Mr_final - realizable moments after initial flux computation
%   flux - struct with fields 'x' and 'y' containing flux moments
%   bounds - struct with eigenvalue bounds for different orders and HLL

%% Step 1: Initial flux computation and realizability
[Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);

%% Step 2: Eigenvalue computation for hyperbolicity
[v6xmin, v6xmax, Mr] = eigenvalues6_hyperbolic_3D(Mr, 'x', flag2D, Ma);
[v6ymin, v6ymax, Mr] = eigenvalues6_hyperbolic_3D(Mr, 'y', flag2D, Ma);

%% Step 3: Eigenvalues for HLL
[v5xmin, v5xmax, v5ymin, v5ymax] = hll_bounds_from_M(Mr, idx);

%% Step 4: Package outputs into structured format
flux = struct('x', Mx, 'y', My);

bounds = struct( ...
    'x', struct('v6min', v6xmin, 'v6max', v6xmax, 'v5min', v5xmin, 'v5max', v5xmax), ...
    'y', struct('v6min', v6ymin, 'v6max', v6ymax, 'v5min', v5ymin, 'v5max', v5ymax), ...
    'hll', struct('xmin', min(v5xmin, v6xmin), 'xmax', max(v5xmax, v6xmax), ...
                  'ymin', min(v5ymin, v6ymin), 'ymax', max(v5ymax, v6ymax)) ...
);

Mr_final = Mr;

end