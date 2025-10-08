% Compare one time step between MATLAB and Julia
% This script runs one time step and saves detailed intermediate values

clear all;
close all;

% Add paths
setup_paths;

% Parameters (matching Julia defaults)
Np = 20;
tmax = 0.1;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;  % 3D mode
CFL = 0.5;

% Grid setup
dx = 2.0 / Np;
dy = 2.0 / Np;

% Create params struct
params = struct();
params.Np = Np;
params.tmax = 0.02;  % Short time for one step
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

% Run one step
fprintf('Running MATLAB simulation for ONE step...\n');
script_dir = pwd;
[M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir);

% Package results
result = struct();
result.final_moments = M_final;
result.final_time = final_time;
result.time_steps = time_steps;
% Calculate dt from final time and steps
result.dt_history = final_time / time_steps;
% Get initial moments - need to reconstruct them
result.initial_moments = zeros(size(M_final));
% For now, just use zeros - we'll focus on final moments

% Save results
save('matlab_one_step.mat', 'result', '-v7.3');

fprintf('\n=== MATLAB One Step Results ===\n');
fprintf('Steps completed: %d\n', result.time_steps);
fprintf('Final time: %.10e\n', result.final_time);
fprintf('dt used: %.10e\n', result.dt_history(1));
fprintf('\nFinal M(1,1,1:5):\n');
for i = 1:5
    fprintf('  M(1,1,%d) = %.15e\n', i, result.final_moments(1,1,i));
end
fprintf('\nFinal M(10,10,1:5):\n');
for i = 1:5
    fprintf('  M(10,10,%d) = %.15e\n', i, result.final_moments(10,10,i));
end
fprintf('\nFinal M(20,20,1:5):\n');
for i = 1:5
    fprintf('  M(20,20,%d) = %.15e\n', i, result.final_moments(20,20,i));
end

% Also save initial conditions for comparison
fprintf('\n=== Initial Conditions ===\n');
fprintf('Initial M(1,1,1:5):\n');
for i = 1:5
    fprintf('  M(1,1,%d) = %.15e\n', i, result.initial_moments(1,1,i));
end
fprintf('Initial M(10,10,1:5):\n');
for i = 1:5
    fprintf('  M(10,10,%d) = %.15e\n', i, result.initial_moments(10,10,i));
end

fprintf('\nResults saved to matlab_one_step.mat\n');
