% Save detailed intermediate values during step 2
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
params.tmax = 0.021;  % Run step 1 first
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

fprintf('Running MATLAB step 1...\n');
script_dir = pwd;
[M_step1, t1, steps1, grid_out] = simulation_runner(params, script_dir);

fprintf('\n=== After Step 1 ===\n');
fprintf('Time: %.15e, Steps: %d\n', t1, steps1);
fprintf('Sample cell (10,10): M(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
        M_step1(10,10,1), M_step1(10,10,2), M_step1(10,10,3), ...
        M_step1(10,10,4), M_step1(10,10,5));

% Now run step 2
params.tmax = 0.04;
fprintf('\nRunning MATLAB step 2...\n');
[M_step2, t2, steps2, grid_out] = simulation_runner(params, script_dir);

fprintf('\n=== After Step 2 ===\n');
fprintf('Time: %.15e, Steps: %d\n', t2, steps2);
fprintf('Sample cell (10,10): M(1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
        M_step2(10,10,1), M_step2(10,10,2), M_step2(10,10,3), ...
        M_step2(10,10,4), M_step2(10,10,5));

% Save both
result = struct();
result.M_step1 = M_step1;
result.M_step2 = M_step2;
result.t1 = t1;
result.t2 = t2;
result.steps1 = steps1;
result.steps2 = steps2;

save('matlab_step2_detailed.mat', 'result', '-v7.3');
fprintf('\nSaved to matlab_step2_detailed.mat\n');
