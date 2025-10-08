% Compare cell (7,13) after step 4
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

% Create params struct - run to t=0.074 (just after step 4)
params = struct();
params.Np = Np;
params.tmax = 0.074;
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

fprintf('Running MATLAB to t=0.074...\n');
script_dir = pwd;
[M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir);

fprintf('\n=== MATLAB Cell (7,13) after step %d ===\n', time_steps);
fprintf('Time: %.15e\n', final_time);
fprintf('M(7,13,1:10):\n');
for k = 1:10
    fprintf('  M(%d) = %.15e\n', k, M_final(7,13,k));
end

% Save
result = struct();
result.M = M_final;
result.t = final_time;
result.steps = time_steps;
save('matlab_cell713.mat', 'result', '-v7.3');
fprintf('\nSaved to matlab_cell713.mat\n');
