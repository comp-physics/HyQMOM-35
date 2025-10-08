% Force exactly TWO time steps for comparison
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

% Create params struct - set tmax to stop after ~2 steps
params = struct();
params.Np = Np;
params.tmax = 0.04;  % Should give us ~2 steps
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

fprintf('Running MATLAB simulation for exactly TWO steps...\n');
script_dir = pwd;
[M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir);

% Package results
result = struct();
result.final_moments = M_final;
result.final_time = final_time;
result.time_steps = time_steps;

% Save results
save('matlab_exactly_two_steps.mat', 'result', '-v7.3');

fprintf('\n=== MATLAB Results ===\n');
fprintf('Steps completed: %d\n', time_steps);
fprintf('Final time: %.15e\n', final_time);

fprintf('\nFinal M(10,10,1:5):\n');
for i = 1:5
    fprintf('  M(10,10,%d) = %.15e\n', i, M_final(10,10,i));
end

fprintf('\nResults saved\n');
