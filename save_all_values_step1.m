% Save ALL intermediate values after step 1
clear all;
close all;

setup_paths;

% Parameters
Np = 20;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;
CFL = 0.5;
dx = 2.0 / Np;
dy = 2.0 / Np;

params = struct();
params.Np = Np;
params.tmax = 0.021;  % Just past step 1
params.Kn = Kn;
params.Ma = Ma;
params.flag2D = flag2D;
params.CFL = CFL;
params.dx = dx;
params.dy = dy;
params.Nmom = 35;
params.nnmax = 1000;
params.dtmax = 0.025;
params.rhol = 1.0;
params.rhor = 0.01;
params.T = 1.0;
params.r110 = 0.0;
params.r101 = 0.0;
params.r011 = 0.0;
params.symmetry_check_interval = 0;
params.enable_memory_tracking = false;

fprintf('Running MATLAB for step 1...\n');
script_dir = pwd;
[M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir);

% Extract derived quantities at each cell
fprintf('Extracting derived quantities...\n');
[Ny, Nx, Nmom] = size(M_final);

% Allocate arrays for derived quantities
rho = zeros(Ny, Nx);
u = zeros(Ny, Nx);
v = zeros(Ny, Nx);
w = zeros(Ny, Nx);
T_field = zeros(Ny, Nx);
C200_field = zeros(Ny, Nx);
C020_field = zeros(Ny, Nx);
C002_field = zeros(Ny, Nx);

for i = 1:Ny
    for j = 1:Nx
        M = squeeze(M_final(i,j,:));
        
        % Basic quantities
        rho(i,j) = M(1);  % M000
        u(i,j) = M(2) / M(1);  % M100/M000
        v(i,j) = M(6) / M(1);  % M010/M000
        w(i,j) = M(16) / M(1); % M001/M000
        
        % Central moments (variances)
        C200_field(i,j) = M(3) - M(2)^2/M(1);  % M200 - M100^2/M000
        C020_field(i,j) = M(10) - M(6)^2/M(1); % M020 - M010^2/M000
        C002_field(i,j) = M(20) - M(16)^2/M(1); % M002 - M001^2/M000
        
        % Temperature (assuming ideal gas, average of diagonal variances)
        T_field(i,j) = (C200_field(i,j) + C020_field(i,j) + C002_field(i,j)) / 3;
    end
end

% Package everything
result = struct();
result.M = M_final;
result.t = final_time;
result.steps = time_steps;
result.rho = rho;
result.u = u;
result.v = v;
result.w = w;
result.T = T_field;
result.C200 = C200_field;
result.C020 = C020_field;
result.C002 = C002_field;
result.params = params;

save('matlab_step1_all_values.mat', 'result', '-v7.3');

fprintf('\n=== MATLAB Step 1 Complete ===\n');
fprintf('Steps: %d, Time: %.15e\n', time_steps, final_time);
fprintf('\nSample values at cell (10,10):\n');
fprintf('  rho = %.15e\n', rho(10,10));
fprintf('  u   = %.15e\n', u(10,10));
fprintf('  v   = %.15e\n', v(10,10));
fprintf('  w   = %.15e\n', w(10,10));
fprintf('  T   = %.15e\n', T_field(10,10));
fprintf('  M(1:5) = [%.6e, %.6e, %.6e, %.6e, %.6e]\n', ...
        M_final(10,10,1), M_final(10,10,2), M_final(10,10,3), ...
        M_final(10,10,4), M_final(10,10,5));

fprintf('\nSaved to matlab_step1_all_values.mat\n');
