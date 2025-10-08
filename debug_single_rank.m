% Single-rank simulation without MPI to match Julia's 1-rank test
% This directly implements the core time-stepping loop

clear; clc;

% Add paths
addpath('src');
addpath('src/autogen');

% Parameters
Np = 20;
tmax = 0.1;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;
CFL = 0.5;

% Grid
dx = 1.0 / Np;
dy = 1.0 / Np;
xm = (0.5:Np) * dx;
ym = (0.5:Np) * dy;

% Initial conditions
rhol = 1.0;
rhor = 0.01;
T = 1.0;

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

fprintf('Initial M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
    M(1,1,1), M(1,1,2), M(1,1,3), M(1,1,4), M(1,1,5));

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

while t < tmax && nn < nnmax
    nn = nn + 1;
    
    % Compute fluxes and wave speeds
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
    
    % X-direction flux update
    decomp = struct('neighbors', struct('left', -1, 'right', -1, 'up', -1, 'down', -1));
    Mnpx = apply_flux_update(M, Fx, vpxmin, vpxmax, vpxmin, vpxmax, ...
                              Np, Np, 0, dt, dx, decomp, 1);
    
    % Y-direction flux update
    Mnpy = apply_flux_update(M, Fy, vpymin, vpymax, vpymin, vpymax, ...
                              Np, Np, 0, dt, dy, decomp, 2);
    
    % Combine updates (Strang splitting)
    M = Mnpx + Mnpy - M;
    
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
    
    % Apply collision
    for i = 1:Np
        for j = 1:Np
            MOM = squeeze(M(i, j, :));
            MMC = collision35(MOM, dt, Kn);
            Mnp(i, j, :) = MMC;
        end
    end
    M = Mnp;
    
    % Check for NaN/Inf
    if any(isnan(M(:))) || any(isinf(M(:)))
        fprintf('ERROR: NaN/Inf detected at step %d!\n', nn);
        break;
    end
    
    t = t + dt;
    fprintf('Step %3d: t = %.6f, dt = %.6e\n', nn, t, dt);
end

fprintf('\nSimulation completed %d steps\n', nn);
fprintf('Final M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
    M(1,1,1), M(1,1,2), M(1,1,3), M(1,1,4), M(1,1,5));
