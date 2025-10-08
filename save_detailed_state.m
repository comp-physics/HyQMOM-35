% Save detailed state after each step for comparison
clear all;
close all;

setup_paths;

% We need to modify simulation_runner to save intermediate states
% For now, let's create a wrapper that saves after each step

Np = 20;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;
CFL = 0.5;
dx = 2.0 / Np;
dy = 2.0 / Np;

params = struct();
params.Np = Np;
params.tmax = 0.074;  % Run to step 5
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

fprintf('Running MATLAB to save detailed state...\n');
script_dir = pwd;
[M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir);

% Save everything
result = struct();
result.M = M_final;
result.t = final_time;
result.steps = time_steps;
result.params = params;

% Save specific cells for comparison
result.cell_7_13 = squeeze(M_final(7,13,:));
result.cell_10_10 = squeeze(M_final(10,10,:));
result.cell_1_1 = squeeze(M_final(1,1,:));

save('matlab_detailed_state.mat', 'result', '-v7.3');

fprintf('\n=== MATLAB Detailed State ===\n');
fprintf('Steps: %d, Time: %.15e\n', time_steps, final_time);
fprintf('\nCell (7,13) moments:\n');
for k = 1:10
    fprintf('  M[%d] = %.15e\n', k, M_final(7,13,k));
end

fprintf('\nSaved to matlab_detailed_state.mat\n');
