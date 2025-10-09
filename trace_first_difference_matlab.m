% Trace the first ~1e-9 difference in the calculation chain
% Step-by-step logging of ALL intermediate values during first time step

clear all;
close all;

setup_paths;

fprintf('===== TRACING FIRST DIFFERENCE - MATLAB =====\n\n');

% Setup (same as actual simulation)
Np = 20;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;
CFL = 0.5;
dx = 2.0 / Np;
dy = 2.0 / Np;
Nmom = 35;

% Initialize grid with halo
halo = 1;
nx = Np;
ny = Np;
M = zeros(nx+2*halo, ny+2*halo, Nmom);

% Initial conditions (exactly as in simulation)
U0 = 0; V0 = 0; W0 = 0;
T = 1.0;
rhol = 1.0;
rhor = 0.01;
r110 = 0.0;
r101 = 0.0;
r011 = 0.0;

C200 = T;
C020 = T;
C002 = T;
C110 = r110 * sqrt(C200 * C020);
C101 = r101 * sqrt(C200 * C002);
C011 = r011 * sqrt(C020 * C002);

% Background state
Mr_bg = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002);

% Crossing jets states
Uc = Ma / sqrt(2);
Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002);
Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002);

% Jet bounds
Csize = floor(0.1 * Np);
Mint = Np/2 + 1;
Maxt = Np/2 + 1 + Csize;
Minb = Np/2 - Csize;
Maxb = Np/2;

% Fill grid
for ii = 1:nx
    gi = ii;
    for jj = 1:ny
        gj = jj;
        
        Mr = Mr_bg;
        
        if gi >= Minb && gi <= Maxb && gj >= Minb && gj <= Maxb
            Mr = Mb;
        end
        
        if gi >= Mint && gi <= Maxt && gj >= Mint && gj <= Maxt
            Mr = Mt;
        end
        
        M(ii + halo, jj + halo, :) = Mr;
    end
end

fprintf('Initial conditions set\n');
fprintf('Tracking cell (8,9) - one of the cells with differences\n\n');

% Track cell (8,9)
track_i = 8;
track_j = 9;
ih = track_i + halo;
jh = track_j + halo;

fprintf('=== INITIAL STATE ===\n');
M_initial = squeeze(M(ih, jh, :));
fprintf('Cell (%d,%d) moments:\n', track_i, track_j);
fprintf('  M000 = %.15e\n', M_initial(1));
fprintf('  M100 = %.15e\n', M_initial(2));
fprintf('  M200 = %.15e\n', M_initial(3));
fprintf('  M210 = %.15e\n', M_initial(8));
fprintf('  M130 = %.15e\n', M_initial(14));
fprintf('\n');

% Compute first time step dt
fprintf('=== COMPUTING TIME STEP ===\n');

% Compute eigenvalues for all cells (simplified - just track our cell)
vpxmin = zeros(nx, ny);
vpxmax = zeros(nx, ny);
vpymin = zeros(nx, ny);
vpymax = zeros(nx, ny);

for i = 1:nx
    for j = 1:ny
        ih_loop = i + halo;
        jh_loop = j + halo;
        MOM = squeeze(M(ih_loop, jh_loop, :));
        [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
        [vpxmin(i,j), vpxmax(i,j), ~] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
        [vpymin(i,j), vpymax(i,j), ~] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
    end
end

vmax = max([abs(vpxmin(:)); abs(vpxmax(:)); abs(vpymin(:)); abs(vpymax(:))]);
dt = CFL * min(dx, dy) / vmax;

fprintf('vmax = %.15e\n', vmax);
fprintf('dt = %.15e\n', dt);
fprintf('\n');

fprintf('=== FLUX COMPUTATION FOR CELL (%d,%d) ===\n', track_i, track_j);

% Get moments at cell and neighbors for flux computation
M_cell = squeeze(M(ih, jh, :));
M_left = squeeze(M(ih-1, jh, :));
M_right = squeeze(M(ih+1, jh, :));
M_down = squeeze(M(ih, jh-1, :));
M_up = squeeze(M(ih, jh+1, :));

fprintf('Before flux computation:\n');
fprintf('  Center M210 = %.15e\n', M_cell(8));
fprintf('  Left M210   = %.15e\n', M_left(8));
fprintf('  Right M210  = %.15e\n', M_right(8));
fprintf('\n');

% Compute fluxes (x-direction)
[~,~,Fx_cell,~] = Flux_closure35_and_realizable_3D(M_cell, flag2D, Ma);
[~,~,Fx_left,~] = Flux_closure35_and_realizable_3D(M_left, flag2D, Ma);
[~,~,Fx_right,~] = Flux_closure35_and_realizable_3D(M_right, flag2D, Ma);

fprintf('X-direction fluxes:\n');
fprintf('  Fx_left(8)  = %.15e\n', Fx_left(8));
fprintf('  Fx_cell(8)  = %.15e\n', Fx_cell(8));
fprintf('  Fx_right(8) = %.15e\n', Fx_right(8));
fprintf('\n');

% HLL flux at left interface
[Fhll_left, vpl_left, vpr_left] = flux_HLL(M_left, M_cell, vpxmin(track_i-1,track_j), ...
                                             vpxmax(track_i-1,track_j), vpxmin(track_i,track_j), ...
                                             vpxmax(track_i,track_j), Fx_left, Fx_cell, 1);

% HLL flux at right interface  
[Fhll_right, vpl_right, vpr_right] = flux_HLL(M_cell, M_right, vpxmin(track_i,track_j), ...
                                                vpxmax(track_i,track_j), vpxmin(track_i+1,track_j), ...
                                                vpxmax(track_i+1,track_j), Fx_cell, Fx_right, 1);

fprintf('HLL fluxes:\n');
fprintf('  Fhll_left(8)  = %.15e\n', Fhll_left(8));
fprintf('  Fhll_right(8) = %.15e\n', Fhll_right(8));
fprintf('\n');

% Apply flux update
Mnpx = M_cell - dt/dx * (Fhll_right - Fhll_left);

fprintf('After X-flux update:\n');
fprintf('  Mnpx(8) = %.15e\n', Mnpx(8));
fprintf('  Change = %.15e\n', Mnpx(8) - M_cell(8));
fprintf('\n');

% Y-direction (similar)
[~,Fy_cell,~,~] = Flux_closure35_and_realizable_3D(M_cell, flag2D, Ma);
[~,Fy_down,~,~] = Flux_closure35_and_realizable_3D(M_down, flag2D, Ma);
[~,Fy_up,~,~] = Flux_closure35_and_realizable_3D(M_up, flag2D, Ma);

[Fhll_down, ~, ~] = flux_HLL(M_down, M_cell, vpymin(track_i,track_j-1), ...
                               vpymax(track_i,track_j-1), vpymin(track_i,track_j), ...
                               vpymax(track_i,track_j), Fy_down, Fy_cell, 2);

[Fhll_up, ~, ~] = flux_HLL(M_cell, M_up, vpymin(track_i,track_j), ...
                             vpymax(track_i,track_j), vpymin(track_i,track_j+1), ...
                             vpymax(track_i,track_j+1), Fy_cell, Fy_up, 2);

Mnpy = M_cell - dt/dy * (Fhll_up - Fhll_down);

fprintf('After Y-flux update:\n');
fprintf('  Mnpy(8) = %.15e\n', Mnpy(8));
fprintf('  Change = %.15e\n', Mnpy(8) - M_cell(8));
fprintf('\n');

% Combine fluxes
Mnp_combined = Mnpx + Mnpy - M_cell;

fprintf('After combining X and Y fluxes:\n');
fprintf('  Mnp_combined(8) = %.15e\n', Mnp_combined(8));
fprintf('  Change from initial = %.15e\n', Mnp_combined(8) - M_cell(8));
fprintf('\n');

% Apply realizability
fprintf('=== REALIZABILITY CORRECTION ===\n');
[~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mnp_combined, flag2D, Ma);
fprintf('After Flux_closure35_and_realizable_3D:\n');
fprintf('  Mr(8) = %.15e\n', Mr(8));
fprintf('  Change = %.15e\n', Mr(8) - Mnp_combined(8));
fprintf('\n');

[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
fprintf('After eigenvalues6_hyperbolic_3D (x):\n');
fprintf('  Mr(8) = %.15e\n', Mr(8));
fprintf('\n');

[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
fprintf('After eigenvalues6_hyperbolic_3D (y):\n');
fprintf('  Mr(8) = %.15e\n', Mr(8));
fprintf('\n');

[~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
fprintf('After second Flux_closure35_and_realizable_3D:\n');
fprintf('  Mr(8) = %.15e\n', Mr(8));
fprintf('\n');

% Apply collision
fprintf('=== COLLISION OPERATOR ===\n');
MMC = collision35(Mr, dt, Kn);
fprintf('After collision35:\n');
fprintf('  MMC(8) = %.15e\n', MMC(8));
fprintf('  Change = %.15e\n', MMC(8) - Mr(8));
fprintf('\n');

fprintf('=== FINAL STATE ===\n');
fprintf('Cell (%d,%d) after one time step:\n', track_i, track_j);
fprintf('  Initial M210 = %.15e\n', M_cell(8));
fprintf('  Final M210   = %.15e\n', MMC(8));
fprintf('  Total change = %.15e\n', MMC(8) - M_cell(8));
fprintf('\n');

fprintf('===== END MATLAB TRACE =====\n');
