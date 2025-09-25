function [Mr_final, Mx, My, v6xmin, v6xmax, v5xmin, v5xmax, v6ymin, v6ymax, v5ymin, v5ymax, vpxmin, vpxmax, vpymin, vpymax] = process_cell_timestep(MOM, flag2D, Ma, idx)
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
%   Mx, My - flux moments
%   v6xmin, v6xmax, v6ymin, v6ymax - 6th order eigenvalue bounds
%   v5xmin, v5xmax, v5ymin, v5ymax - 5th order eigenvalue bounds  
%   vpxmin, vpxmax, vpymin, vpymax - combined eigenvalue bounds for HLL

%% Step 1: Initial flux computation and realizability
[Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);

%% Step 2: Eigenvalue computation for hyperbolicity
[v6xmin, v6xmax, Mr] = eigenvalues6_hyperbolic_3D(Mr, 'x', flag2D, Ma);
[v6ymin, v6ymax, Mr] = eigenvalues6_hyperbolic_3D(Mr, 'y', flag2D, Ma);

%% Step 3: Eigenvalues for HLL
% 1-D hyqmom for m500 eigenvalues in x direction
MOM5 = Mr(idx.x_moments); % m000 m100 m200 m300 m400
[~, v5xmin, v5xmax] = closure_and_eigenvalues(MOM5);

% 1-D hyqmom for m050 eigenvalues in y direction  
MOM5 = Mr(idx.y_moments); % m000 m010 m020 m030 m040
[~, v5ymin, v5ymax] = closure_and_eigenvalues(MOM5);

% Combined eigenvalue bounds for HLL
vpxmin = min(v5xmin, v6xmin);
vpxmax = max(v5xmax, v6xmax);
vpymin = min(v5ymin, v6ymin);
vpymax = max(v5ymax, v6ymax);

Mr_final = Mr;

end