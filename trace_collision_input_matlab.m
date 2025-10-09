% Simplified trace: What does the collision operator see?
% This traces through ONE cell to see where nonzero cross-moments appear

clear all;
close all;

setup_paths;

fprintf('===== TRACE COLLISION INPUT - MATLAB =====\n\n');

% Run one timestep with instrumented collision operator
Np = 20;
Kn = 1.0;
Ma = 0.0;
tmax = 0.1;  % Just one step
flag2D = 0;
CFL = 0.5;

% Track specific cell
global TRACK_CELL;
TRACK_CELL = struct('i', 8, 'j', 9, 'log', []);

params = struct('Np', Np, 'Kn', Kn, 'Ma', Ma, 'tmax', tmax, 'flag2D', flag2D, 'CFL', CFL);

% Patch collision35 to log inputs
original_collision35 = @collision35;
collision35_instrumented = @(M, dt, Kn) collision35_logged(M, dt, Kn, original_collision35);

% Run simulation (we'll manually instrument it)
fprintf('Running simulation for one step...\n');

% Initialize
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
            Mr = Mt;
        end
        M(ii + halo, jj + halo, :) = Mr;
    end
end

track_i = 8;
track_j = 9;
ih = track_i + halo;
jh = track_j + halo;

fprintf('Cell (%d,%d) initial state:\n', track_i, track_j);
fprintf('  M210 = %.15e\n', M(ih, jh, 8));
fprintf('  M130 = %.15e\n\n', M(ih, jh, 14));

% Compute dt
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
fprintf('dt = %.15e\n\n', dt);

% Flux computation
Mx = zeros(nx+2*halo, ny+2*halo, Nmom);
My = zeros(nx+2*halo, ny+2*halo, Nmom);
Fy = zeros(nx+2*halo, ny+2*halo, Nmom);
Fx = zeros(nx+2*halo, ny+2*halo, Nmom);

for i = 1:nx
    for j = 1:ny
        ih_loop = i + halo;
        jh_loop = j + halo;
        MOM = squeeze(M(ih_loop, jh_loop, :));
        [~,Fy(ih_loop,jh_loop,:),Fx(ih_loop,jh_loop,:),Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
        Mx(ih_loop, jh_loop, :) = Mr;
        My(ih_loop, jh_loop, :) = Mr;
    end
end

% Halo exchange (periodic)
M(:,1,:) = M(:,ny+1,:);
M(:,ny+2*halo,:) = M(:,halo+1,:);
M(1,:,:) = M(nx+1,:,:);
M(nx+2*halo,:,:) = M(halo+1,:,:);

% X-direction update
Mnpx = pas_HLL(squeeze(Mx(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               squeeze(Fx(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               dt, dx, vpxmin, vpxmax, true, true);

% Y-direction update
Mnpy = pas_HLL(squeeze(My(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               squeeze(Fy(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               dt, dy, vpymin, vpymax, true, true);

% Combine
Mnp_temp = zeros(nx+2*halo, ny+2*halo, Nmom);
for i = 1:nx
    for j = 1:ny
        Mnp_temp(i+halo, j+halo, :) = Mnpx(i,j,:) + Mnpy(i,j,:) - squeeze(M(i+halo, j+halo, :));
    end
end

fprintf('After flux update at cell (%d,%d):\n', track_i, track_j);
fprintf('  M210 = %.15e\n', Mnp_temp(ih, jh, 8));
fprintf('  M130 = %.15e\n\n', Mnp_temp(ih, jh, 14));

% Realizability
Mnp = zeros(nx+2*halo, ny+2*halo, Nmom);
v6xmin = zeros(nx, ny);
v6xmax = zeros(nx, ny);
v6ymin = zeros(nx, ny);
v6ymax = zeros(nx, ny);

for i = 1:nx
    for j = 1:ny
        ih_loop = i + halo;
        jh_loop = j + halo;
        MOM = squeeze(Mnp_temp(ih_loop, jh_loop, :));
        
        [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
        [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
        [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
        [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
        
        Mnp(ih_loop, jh_loop, :) = Mr;
    end
end

fprintf('After realizability at cell (%d,%d):\n', track_i, track_j);
fprintf('  M210 = %.15e\n', Mnp(ih, jh, 8));
fprintf('  M130 = %.15e\n\n', Mnp(ih, jh, 14));

% Halo exchange
Mnp(:,1,:) = Mnp(:,ny+1,:);
Mnp(:,ny+2*halo,:) = Mnp(:,halo+1,:);
Mnp(1,:,:) = Mnp(nx+1,:,:);
Mnp(nx+2*halo,:,:) = Mnp(halo+1,:,:);

% Collision
fprintf('=== BEFORE COLLISION ===\n');
M_before_collision = squeeze(Mnp(ih, jh, :));
fprintf('  M210 = %.15e\n', M_before_collision(8));
fprintf('  M130 = %.15e\n\n', M_before_collision(14));

MMC = collision35(M_before_collision, dt, Kn);

fprintf('=== AFTER COLLISION ===\n');
fprintf('  M210 = %.15e\n', MMC(8));
fprintf('  M130 = %.15e\n\n', MMC(14));

fprintf('=== FINAL ANSWER ===\n');
fprintf('Cell (%d,%d) moments going into collision operator:\n', track_i, track_j);
fprintf('  M210 = %.15e\n', M_before_collision(8));
fprintf('  M130 = %.15e\n', M_before_collision(14));

fprintf('\n===== END MATLAB TRACE =====\n');
