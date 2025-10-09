% Save complete state at end of step 4 for Julia comparison
clear all;
close all;

setup_paths;

fprintf('===== SAVING STEP 4 STATE FOR COMPARISON =====\n\n');

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
halo = 2;  % Matching Julia's halo size
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

fprintf('Running simulation to step 4...\n');

% Run exactly 4 steps
t = 0.0;
for step = 1:4
    fprintf('  Step %d...', step);
    
    % Compute wave speeds and fluxes (simplified - just get dt)
    vpxmin = zeros(nx, ny);
    vpxmax = zeros(nx, ny);
    vpymin = zeros(nx, ny);
    vpymax = zeros(nx, ny);
    
    Fx = zeros(nx+2*halo, ny+2*halo, Nmom);
    Fy = zeros(nx+2*halo, ny+2*halo, Nmom);
    Mnp = zeros(nx+2*halo, ny+2*halo, Nmom);
    
    % Compute fluxes for interior cells
    for i = 1:nx
        for j = 1:ny
            ih = i + halo;
            jh = j + halo;
            MOM = squeeze(M(ih, jh, :));
            
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [vpxmin(i,j), vpxmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [vpymin(i,j), vpymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            
            Fx(ih, jh, :) = Mx;
            Fy(ih, jh, :) = My;
            Mnp(ih, jh, :) = Mr;
        end
    end
    
    M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
    
    % Apply boundary conditions (copy)
    M(1:halo, :, :) = M(halo+1:halo+halo, :, :);
    M(halo+nx+1:end, :, :) = M(halo+nx-halo+1:halo+nx, :, :);
    M(:, 1:halo, :) = M(:, halo+1:halo+halo, :);
    M(:, halo+ny+1:end, :) = M(:, halo+ny-halo+1:halo+ny, :);
    
    % Also exchange Fx, Fy boundaries
    Fx(1:halo, :, :) = Fx(halo+1:halo+halo, :, :);
    Fx(halo+nx+1:end, :, :) = Fx(halo+nx-halo+1:halo+nx, :, :);
    Fx(:, 1:halo, :) = Fx(:, halo+1:halo+halo, :);
    Fx(:, halo+ny+1:end, :) = Fx(:, halo+ny-halo+1:halo+ny, :);
    
    Fy(1:halo, :, :) = Fy(halo+1:halo+halo, :, :);
    Fy(halo+nx+1:end, :, :) = Fy(halo+nx-halo+1:halo+nx, :, :);
    Fy(:, 1:halo, :) = Fy(:, halo+1:halo+halo, :);
    Fy(:, halo+ny+1:end, :) = Fy(:, halo+ny-halo+1:halo+ny, :);
    
    % Compute time step
    vmax = max([abs(vpxmin(:)); abs(vpxmax(:)); abs(vpymin(:)); abs(vpymax(:))]);
    dt = CFL * min(dx, dy) / vmax;
    t = t + dt;
    
    % Apply flux updates (simplified HLL)
    for i = 1:nx
        for j = 1:ny
            ih = i + halo;
            jh = j + halo;
            
            % X-direction
            Mnpx = squeeze(M(ih, jh, :)) - dt/dx * (squeeze(Fx(ih+1, jh, :)) - squeeze(Fx(ih, jh, :)));
            
            % Y-direction
            Mnpy = squeeze(M(ih, jh, :)) - dt/dy * (squeeze(Fy(ih, jh+1, :)) - squeeze(Fy(ih, jh, :)));
            
            % Combine
            M_combined = Mnpx + Mnpy - squeeze(M(ih, jh, :));
            
            Mnp(ih, jh, :) = reshape(M_combined, [1,1,Nmom]);
        end
    end
    
    M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
    
    % Apply boundary conditions again
    M(1:halo, :, :) = M(halo+1:halo+halo, :, :);
    M(halo+nx+1:end, :, :) = M(halo+nx-halo+1:halo+nx, :, :);
    M(:, 1:halo, :) = M(:, halo+1:halo+halo, :);
    M(:, halo+ny+1:end, :) = M(:, halo+ny-halo+1:halo+ny, :);
    
    % Realizability
    for i = 1:nx
        for j = 1:ny
            ih = i + halo;
            jh = j + halo;
            MOM = squeeze(M(ih, jh, :));
            
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            
            Mnp(ih, jh, :) = Mr;
        end
    end
    
    M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
    
    % Collision
    for i = 1:nx
        for j = 1:ny
            ih = i + halo;
            jh = j + halo;
            MOM = squeeze(M(ih, jh, :));
            MMC = collision35(MOM, dt, Kn);
            Mnp(ih, jh, :) = MMC;
        end
    end
    
    M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
    
    % Apply boundary conditions final
    M(1:halo, :, :) = M(halo+1:halo+halo, :, :);
    M(halo+nx+1:end, :, :) = M(halo+nx-halo+1:halo+nx, :, :);
    M(:, 1:halo, :) = M(:, halo+1:halo+halo, :);
    M(:, halo+ny+1:end, :) = M(:, halo+ny-halo+1:halo+ny, :);
    
    fprintf(' t=%.6f, dt=%.6e\n', t, dt);
end

fprintf('\nStep 4 complete. Saving state...\n');

% Save everything needed for step 5
save('matlab_step4_complete.mat', 'M', 'Mnp', 'Fx', 'Fy', ...
     'vpxmin', 'vpxmax', 'vpymin', 'vpymax', ...
     't', 'dt', 'nx', 'ny', 'halo', 'Nmom', ...
     'Np', 'Kn', 'Ma', 'flag2D', 'CFL', 'dx', 'dy', '-v7.3');

% Also save cell (7,13) specifically for detailed tracking
ih_track = 7 + halo;
jh_track = 13 + halo;
M_cell_713 = squeeze(M(ih_track, jh_track, :));

fprintf('\nCell (7,13) at end of step 4:\n');
fprintf('  M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n', M_cell_713(1:5));
fprintf('  M[8] = %.6e (M210)\n', M_cell_713(8));
fprintf('  M[14] = %.6e (M130)\n', M_cell_713(14));

save('matlab_step4_cell713.mat', 'M_cell_713', 'ih_track', 'jh_track');

fprintf('\n✅ Saved to matlab_step4_complete.mat\n');
fprintf('✅ Saved cell (7,13) to matlab_step4_cell713.mat\n');

