% Trace ALL intermediate values during first time step
% Find exactly where MATLAB and Julia first diverge

clear all;
close all;

setup_paths;

fprintf('===== TRACING INTERMEDIATE VALUES - MATLAB =====\n\n');

% Setup
Np = 20;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;
CFL = 0.5;
dx = 2.0 / Np;
dy = 2.0 / Np;
Nmom = 35;

% Initialize
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

fprintf('=== INITIAL CONDITION PARAMETERS ===\n');
fprintf('C200 = %.15e\n', C200);
fprintf('C110 = %.15e\n', C110);
fprintf('C101 = %.15e\n', C101);
fprintf('C020 = %.15e\n', C020);
fprintf('C011 = %.15e\n', C011);
fprintf('C002 = %.15e\n\n', C002);

Mr_bg = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002);
Uc = Ma / sqrt(2);
Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002);
Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002);

fprintf('=== INITIAL MOMENTS ===\n');
fprintf('Mr_bg(8)  = %.15e\n', Mr_bg(8));
fprintf('Mr_bg(14) = %.15e\n', Mr_bg(14));
fprintf('Mt(8)  = %.15e\n', Mt(8));
fprintf('Mt(14) = %.15e\n', Mt(14));
fprintf('Mb(8)  = %.15e\n', Mb(8));
fprintf('Mb(14) = %.15e\n\n', Mb(14));

% Fill grid
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

% Track cell (8,9)
track_i = 8;
track_j = 9;
ih = track_i + halo;
jh = track_j + halo;

fprintf('=== CELL (%d,%d) INITIAL STATE ===\n', track_i, track_j);
M_init = squeeze(M(ih, jh, :));
fprintf('M(8)  = %.15e\n', M_init(8));
fprintf('M(14) = %.15e\n\n', M_init(14));

% Compute dt
vpxmin = zeros(nx, ny);
vpxmax = zeros(nx, ny);
vpymin = zeros(nx, ny);
vpymax = zeros(nx, ny);

fprintf('=== COMPUTING EIGENVALUES (first few cells) ===\n');
for i = 1:min(3, nx)
    for j = 1:min(3, ny)
        ih_loop = i + halo;
        jh_loop = j + halo;
        MOM = squeeze(M(ih_loop, jh_loop, :));
        
        [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
        [vpxmin(i,j), vpxmax(i,j), ~] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
        [vpymin(i,j), vpymax(i,j), ~] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
        
        if i == 1 && j == 1
            fprintf('Cell (1,1):\n');
            fprintf('  vpxmin = %.15e, vpxmax = %.15e\n', vpxmin(i,j), vpxmax(i,j));
            fprintf('  vpymin = %.15e, vpymax = %.15e\n', vpymin(i,j), vpymax(i,j));
        end
    end
end

% Complete eigenvalue computation for all cells
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

fprintf('\nvmax = %.15e\n', vmax);
fprintf('dt   = %.15e\n\n', dt);

% Flux computation
fprintf('=== FLUX COMPUTATION (cell %d,%d) ===\n', track_i, track_j);
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
        
        if i == track_i && j == track_j
            fprintf('After Flux_closure35:\n');
            fprintf('  Mr(8)  = %.15e\n', Mr(8));
            fprintf('  Mr(14) = %.15e\n', Mr(14));
            fprintf('  Fx(8)  = %.15e\n', Fx(ih_loop,jh_loop,8));
            fprintf('  Fy(8)  = %.15e\n\n', Fy(ih_loop,jh_loop,8));
        end
    end
end

% Halo exchange
M(:,1,:) = M(:,ny+1,:);
M(:,ny+2*halo,:) = M(:,halo+1,:);
M(1,:,:) = M(nx+1,:,:);
M(nx+2*halo,:,:) = M(halo+1,:,:);

% Flux updates
fprintf('=== FLUX UPDATES ===\n');
Mnpx = pas_HLL(squeeze(Mx(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               squeeze(Fx(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               dt, dx, vpxmin, vpxmax, true, true);

Mnpy = pas_HLL(squeeze(My(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               squeeze(Fy(halo+1:halo+nx, halo+1:halo+ny, :)), ...
               dt, dy, vpymin, vpymax, true, true);

fprintf('Cell (%d,%d) after flux updates:\n', track_i, track_j);
fprintf('  Mnpx(8) = %.15e\n', Mnpx(track_i, track_j, 8));
fprintf('  Mnpy(8) = %.15e\n\n', Mnpy(track_i, track_j, 8));

% Combine
Mnp_temp = zeros(nx+2*halo, ny+2*halo, Nmom);
for i = 1:nx
    for j = 1:ny
        Mnp_temp(i+halo, j+halo, :) = Mnpx(i,j,:) + Mnpy(i,j,:) - squeeze(M(i+halo, j+halo, :));
    end
end

fprintf('After combining fluxes:\n');
fprintf('  Mnp_temp(%d,%d,8) = %.15e\n\n', ih, jh, Mnp_temp(ih, jh, 8));

% Realizability
fprintf('=== REALIZABILITY CORRECTIONS ===\n');
Mnp = zeros(nx+2*halo, ny+2*halo, Nmom);

for i = 1:nx
    for j = 1:ny
        ih_loop = i + halo;
        jh_loop = j + halo;
        MOM = squeeze(Mnp_temp(ih_loop, jh_loop, :));
        
        [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
        [~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
        [~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
        [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
        
        Mnp(ih_loop, jh_loop, :) = Mr;
        
        if i == track_i && j == track_j
            fprintf('Cell (%d,%d) after realizability:\n', i, j);
            fprintf('  Mr(8)  = %.15e\n', Mr(8));
            fprintf('  Mr(14) = %.15e\n\n', Mr(14));
        end
    end
end

% Halo exchange
Mnp(:,1,:) = Mnp(:,ny+1,:);
Mnp(:,ny+2*halo,:) = Mnp(:,halo+1,:);
Mnp(1,:,:) = Mnp(nx+1,:,:);
Mnp(nx+2*halo,:,:) = Mnp(halo+1,:,:);

% Bulk assignment
M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);

% Collision
fprintf('=== COLLISION OPERATOR ===\n');
M_before_coll = squeeze(M(ih, jh, :));
fprintf('Before collision:\n');
fprintf('  M(8)  = %.15e\n', M_before_coll(8));
fprintf('  M(14) = %.15e\n', M_before_coll(14));

% Extract covariance inputs to collision
rho_coll = M_before_coll(1);
umean_coll = M_before_coll(2) / rho_coll;
vmean_coll = M_before_coll(6) / rho_coll;
wmean_coll = M_before_coll(16) / rho_coll;
C200_coll = M_before_coll(3)/rho_coll - umean_coll^2;
C020_coll = M_before_coll(10)/rho_coll - vmean_coll^2;
C002_coll = M_before_coll(20)/rho_coll - wmean_coll^2;
Theta_coll = (C200_coll + C020_coll + C002_coll) / 3;

fprintf('\nCovariance inputs to collision:\n');
fprintf('  rho   = %.15e\n', rho_coll);
fprintf('  umean = %.15e\n', umean_coll);
fprintf('  C200  = %.15e\n', C200_coll);
fprintf('  C020  = %.15e\n', C020_coll);
fprintf('  C002  = %.15e\n', C002_coll);
fprintf('  Theta = %.15e\n\n', Theta_coll);

MMC = collision35(M_before_coll, dt, Kn);

fprintf('After collision:\n');
fprintf('  M(8)  = %.15e\n', MMC(8));
fprintf('  M(14) = %.15e\n\n', MMC(14));

% Save all intermediate values
save('matlab_intermediate_trace.mat', 'M_init', 'M_before_coll', 'MMC', ...
     'dt', 'vmax', 'rho_coll', 'umean_coll', 'C200_coll', 'C020_coll', ...
     'C002_coll', 'Theta_coll', 'track_i', 'track_j');

fprintf('Saved to matlab_intermediate_trace.mat\n');
fprintf('\n===== END MATLAB TRACE =====\n');
