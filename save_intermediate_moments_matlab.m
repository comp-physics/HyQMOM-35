% Save intermediate moments during first timestep
% Focus on cell (8,9) where differences occur

clear all;
close all;

setup_paths;

fprintf('===== SAVING INTERMEDIATE MOMENTS - MATLAB =====\n\n');

% Run full simulation and capture M after each stage
params = struct('Np', 20, 'Kn', 1.0, 'Ma', 0.0, 'tmax', 0.05, 'flag2D', 0, 'CFL', 0.5);

fprintf('Running simulation for one timestep...\n');
result = simulation_runner(params);

% The result contains final moments, but we need intermediate stages
% Let's manually run one step with logging

% Initialize same as simulation
Np = 20;
flag2D = 0;
Ma = 0.0;
Kn = 1.0;
CFL = 0.5;

% Grid setup
dx = 2.0 / Np;
dy = 2.0 / Np;
Nmom = 35;
halo = 1;
nx = Np;
ny = Np;
M = zeros(nx+2*halo, ny+2*halo, Nmom);

% Initial conditions
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

Mr_bg = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002);
Uc = Ma / sqrt(2);
Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002);
Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002);

Csize = floor(0.1 * Np);
Mint = Np/2 + 1;
Maxt = Np/2 + 1 + Csize;
Minb = Np/2 - Csize;
Maxb = Np/2;

for ii = 1:nx
    for jj = 1:ny
        Mr = Mr_bg;
        if ii >= Minb && ii <= Maxb && jj >= Minb && jj <= Maxb
            Mr = Mb;
        end
        if ii >= Mint && ii <= Maxt && jj >= Mint && jj <= Maxt
            Mr = Mt
;
        end
        M(ii + halo, jj + halo, :) = Mr;
    end
end

track_i = 8;
track_j = 9;
ih = track_i + halo;
jh = track_j + halo;

stages = struct();
stages.initial = squeeze(M(ih, jh, :));

fprintf('Stage: Initial\n');
fprintf('  M210 = %.15e, M130 = %.15e\n\n', stages.initial(8), stages.initial(14));

% Compute dt (same as simulation)
[~, ~, dt] = compute_timestep_and_eigenvalues(M, dx, dy, CFL, flag2D, Ma);

fprintf('dt = %.15e\n\n', dt);

% Apply one timestep manually
[M_after_step, ~] = apply_one_timestep(M, dx, dy, dt, Kn, flag2D, Ma);

stages.after_step = squeeze(M_after_step(ih, jh, :));

fprintf('Stage: After full timestep\n');
fprintf('  M210 = %.15e, M130 = %.15e\n\n', stages.after_step(8), stages.after_step(14));

% Save to file
save('matlab_intermediate_moments.mat', 'stages', 'track_i', 'track_j', 'dt');

fprintf('Saved to matlab_intermediate_moments.mat\n');
fprintf('===== END MATLAB =====\n');

% Helper functions
function [M_new, vmax] = apply_one_timestep(M, dx, dy, dt, Kn, flag2D, Ma)
    [~, ~, ~, M_new, vmax] = timestep_explicit(M, dx, dy, dt, Kn, flag2D, Ma);
end

function [vmax, eigenvals, dt] = compute_timestep_and_eigenvalues(M, dx, dy, CFL, flag2D, Ma)
    halo = 1;
    nx = size(M, 1) - 2*halo;
    ny = size(M, 2) - 2*halo;
    
    vpxmin = zeros(nx, ny);
    vpxmax = zeros(nx, ny);
    vpymin = zeros(nx, ny);
    vpymax = zeros(nx, ny);
    
    for i = 1:nx
        for j = 1:ny
            ih = i + halo;
            jh = j + halo;
            MOM = squeeze(M(ih, jh, :));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [vpxmin(i,j), vpxmax(i,j), ~] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [vpymin(i,j), vpymax(i,j), ~] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
        end
    end
    
    vmax = max([abs(vpxmin(:)); abs(vpxmax(:)); abs(vpymin(:)); abs(vpymax(:))]);
    dt = CFL * min(dx, dy) / vmax;
    eigenvals = struct('vpxmin', vpxmin, 'vpxmax', vpxmax, 'vpymin', vpymin, 'vpymax', vpymax);
end
