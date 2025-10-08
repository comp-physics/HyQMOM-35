% Debug script to trace exactly what happens in MATLAB simulation
% Run this to compare with Julia step-by-step

clear; clc;

% Parameters matching the golden file
Np = 20;
tmax = 0.1;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;
CFL = 0.5;

% Initial condition parameters (crossing jets)
rhol = 1.0;
rhor = 0.01;  % Low density background
T = 1.0;

% Grid
dx = 1.0 / Np;
dy = 1.0 / Np;
xm = (0.5:Np) * dx;
ym = (0.5:Np) * dy;

% Initialize moments
M = zeros(Np, Np, 35);
for i = 1:Np
    for j = 1:Np
        x = xm(i);
        y = ym(j);
        
        % Crossing jets IC
        if (x >= 0.25 && x <= 0.75) && (y < 0.25 || y > 0.75)
            rho = rhol;
        elseif (y >= 0.25 && y <= 0.75) && (x < 0.25 || x > 0.75)
            rho = rhol;
        else
            rho = rhor;
        end
        
        M_cell = InitializeM4_35(rho, 0, 0, 0, T, 0, 0, T, 0, T);
        M(i, j, :) = M_cell;
    end
end

fprintf('=== INITIAL STATE ===\n');
fprintf('M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', M(1,1,1), M(1,1,2), M(1,1,3), M(1,1,4), M(1,1,5));
fprintf('M(10,10,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', M(10,10,1), M(10,10,2), M(10,10,3), M(10,10,4), M(10,10,5));
fprintf('\n');

% Time stepping
t = 0;
nn = 0;
nnmax = 1000;
dtmax = CFL * min(dx, dy);

% Preallocate
Fx = zeros(Np, Np, 35);
Fy = zeros(Np, Np, 35);
Mnp = zeros(Np, Np, 35);
v6xmin = zeros(Np, Np);
v6xmax = zeros(Np, Np);
v6ymin = zeros(Np, Np);
v6ymax = zeros(Np, Np);
v5xmin = zeros(Np, Np);
v5xmax = zeros(Np, Np);
v5ymin = zeros(Np, Np);
v5ymax = zeros(Np, Np);
vpxmin = zeros(Np, Np);
vpxmax = zeros(Np, Np);
vpymin = zeros(Np, Np);
vpymax = zeros(Np, Np);

% Run just a few steps with detailed output
max_debug_steps = 7;

while t < tmax && nn < nnmax && nn < max_debug_steps
    nn = nn + 1;
    
    fprintf('=== STEP %d START (t=%.15e) ===\n', nn, t);
    
    % Compute fluxes for cell (1,1) with detailed output
    i = 1; j = 1;
    MOM = squeeze(M(i, j, :));
    
    fprintf('Cell (1,1) before flux computation:\n');
    fprintf('  M(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', MOM(1), MOM(2), MOM(3), MOM(4), MOM(5));
    
    % First realizability call
    [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
    fprintf('  After Flux_closure35 #1:\n');
    fprintf('    Mr(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', Mr(1), Mr(2), Mr(3), Mr(4), Mr(5));
    
    % Eigenvalues X
    [v6xmin_val, v6xmax_val, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
    fprintf('  After eigenvalues6 X:\n');
    fprintf('    v6xmin=%.15e, v6xmax=%.15e\n', v6xmin_val, v6xmax_val);
    fprintf('    Mr(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', Mr(1), Mr(2), Mr(3), Mr(4), Mr(5));
    
    % Eigenvalues Y
    [v6ymin_val, v6ymax_val, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
    fprintf('  After eigenvalues6 Y:\n');
    fprintf('    v6ymin=%.15e, v6ymax=%.15e\n', v6ymin_val, v6ymax_val);
    fprintf('    Mr(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', Mr(1), Mr(2), Mr(3), Mr(4), Mr(5));
    
    % Second realizability call
    [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
    fprintf('  After Flux_closure35 #2:\n');
    fprintf('    Mx(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', Mx(1), Mx(2), Mx(3), Mx(4), Mx(5));
    fprintf('    My(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', My(1), My(2), My(3), My(4), My(5));
    fprintf('    Mr(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', Mr(1), Mr(2), Mr(3), Mr(4), Mr(5));
    
    % Closure and eigenvalues
    [~, v5xmin_val, v5xmax_val] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
    [~, v5ymin_val, v5ymax_val] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
    fprintf('  Closure eigenvalues:\n');
    fprintf('    v5xmin=%.15e, v5xmax=%.15e\n', v5xmin_val, v5xmax_val);
    fprintf('    v5ymin=%.15e, v5ymax=%.15e\n', v5ymin_val, v5ymax_val);
    
    % Check for NaN/Inf
    if any(isnan(Mr)) || any(isinf(Mr))
        fprintf('  WARNING: Mr contains NaN or Inf!\n');
    end
    
    % Now compute for all cells (abbreviated output)
    for i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i, j, :));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            
            Fx(i, j, :) = Mx;
            Fy(i, j, :) = My;
            Mnp(i, j, :) = Mr;
            
            [~, v5xmin(i,j), v5xmax(i,j)] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
            [~, v5ymin(i,j), v5ymax(i,j)] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
            
            vpxmin(i,j) = min(v5xmin(i,j), v6xmin(i,j));
            vpxmax(i,j) = max(v5xmax(i,j), v6xmax(i,j));
            vpymin(i,j) = min(v5ymin(i,j), v6ymin(i,j));
            vpymax(i,j) = max(v5ymax(i,j), v6ymax(i,j));
        end
    end
    M = Mnp;
    
    % Compute dt
    vmax = max([abs(vpxmax(:)); abs(vpxmin(:)); abs(vpymax(:)); abs(vpymin(:))]);
    dt = min(CFL*dx/vmax, dtmax);
    dt = min(dt, tmax-t);
    
    fprintf('  vmax=%.15e, dt=%.15e\n', vmax, dt);
    
    % Apply flux update (simplified for 1 rank, no boundaries)
    % X-direction
    for j = 1:Np
        for i = 1:Np
            if i == 1
                % Left boundary
                Mnpx(i,j,:) = M(i,j,:);
            elseif i == Np
                % Right boundary  
                Mnpx(i,j,:) = M(i,j,:);
            else
                % Interior
                W_L = squeeze(M(i-1,j,:));
                W_R = squeeze(M(i,j,:));
                F_L = squeeze(Fx(i-1,j,:));
                F_R = squeeze(Fx(i,j,:));
                vp_L = vpxmax(i-1,j);
                vp_R = vpxmin(i,j);
                
                [flux, ~] = flux_HLL(W_L, W_R, F_L, F_R, vp_L, vp_R);
                
                W_L = squeeze(M(i,j,:));
                W_R = squeeze(M(i+1,j,:));
                F_L = squeeze(Fx(i,j,:));
                F_R = squeeze(Fx(i+1,j,:));
                vp_L = vpxmax(i,j);
                vp_R = vpxmin(i+1,j);
                
                [flux_p1, ~] = flux_HLL(W_L, W_R, F_L, F_R, vp_L, vp_R);
                
                Mnpx(i,j,:) = M(i,j,:) - (dt/dx) * (flux_p1 - flux);
            end
        end
    end
    
    % Y-direction (similar)
    for i = 1:Np
        for j = 1:Np
            if j == 1 || j == Np
                Mnpy(i,j,:) = M(i,j,:);
            else
                W_L = squeeze(M(i,j-1,:));
                W_R = squeeze(M(i,j,:));
                F_L = squeeze(Fy(i,j-1,:));
                F_R = squeeze(Fy(i,j,:));
                vp_L = vpymax(i,j-1);
                vp_R = vpymin(i,j);
                
                [flux, ~] = flux_HLL(W_L, W_R, F_L, F_R, vp_L, vp_R);
                
                W_L = squeeze(M(i,j,:));
                W_R = squeeze(M(i,j+1,:));
                F_L = squeeze(Fy(i,j,:));
                F_R = squeeze(Fy(i,j+1,:));
                vp_L = vpymax(i,j);
                vp_R = vpymin(i,j+1);
                
                [flux_p1, ~] = flux_HLL(W_L, W_R, F_L, F_R, vp_L, vp_R);
                
                Mnpy(i,j,:) = M(i,j,:) - (dt/dy) * (flux_p1 - flux);
            end
        end
    end
    
    % Combine (Strang splitting)
    M = Mnpx + Mnpy - M;
    
    fprintf('  After flux update, M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
        M(1,1,1), M(1,1,2), M(1,1,3), M(1,1,4), M(1,1,5));
    
    % Enforce realizability
    for i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i, j, :));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            Mnp(i, j, :) = Mr;
        end
    end
    M = Mnp;
    
    fprintf('  After realizability, M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
        M(1,1,1), M(1,1,2), M(1,1,3), M(1,1,4), M(1,1,5));
    
    % Apply collision
    for i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i, j, :));
            MMC = collision35(MOM, dt, Kn);
            Mnp(i, j, :) = MMC;
        end
    end
    M = Mnp;
    
    fprintf('  After collision, M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
        M(1,1,1), M(1,1,2), M(1,1,3), M(1,1,4), M(1,1,5));
    
    % Check for NaN/Inf
    if any(isnan(M(:))) || any(isinf(M(:)))
        fprintf('  ERROR: M contains NaN or Inf after step %d!\n', nn);
        fprintf('  Number of NaN cells: %d\n', sum(any(isnan(M), 3), 'all'));
        fprintf('  Number of Inf cells: %d\n', sum(any(isinf(M), 3), 'all'));
        break;
    end
    
    t = t + dt;
    fprintf('=== STEP %d END (t=%.15e) ===\n\n', nn, t);
end

fprintf('Simulation completed %d steps\n', nn);
fprintf('Final M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
    M(1,1,1), M(1,1,2), M(1,1,3), M(1,1,4), M(1,1,5));
