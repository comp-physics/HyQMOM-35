% Debug realizability at step 4, cell (7,13)
% We need to instrument the MATLAB code to print what happens

clear all;
close all;

setup_paths;

% Run to just before step 4 completes
Np = 20;
Kn = 1.0;
Ma = 0.0;
flag2D = 0;
CFL = 0.5;
dx = 2.0 / Np;
dy = 2.0 / Np;

params = struct();
params.Np = Np;
params.tmax = 0.073;  % Stop just before step 4 completes
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

fprintf('Running MATLAB to just before step 4...\n');
script_dir = pwd;
[M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir);

fprintf('\n=== MATLAB at t=%.6f (step %d) ===\n', final_time, time_steps);
fprintf('Cell (7,13) M[1:5]:\n');
for k = 1:5
    fprintf('  M[%d] = %.15e\n', k, M_final(7,13,k));
end

% Now manually test what realizability does to this cell
MOM = squeeze(M_final(7,13,:));
fprintf('\n=== Testing Flux_closure35_and_realizable_3D ===\n');
fprintf('Input M[3] = %.15e\n', MOM(3));

[Mx, My, Mz, Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);

fprintf('After realizability M[3] = %.15e\n', Mr(3));
fprintf('Change: %.15e (%.2fx)\n', Mr(3) - MOM(3), Mr(3) / MOM(3));

% Save for Julia comparison
result = struct();
result.MOM_before = MOM;
result.Mr_after = Mr;
result.Mx = Mx;
result.My = My;
result.Mz = Mz;
save('matlab_realizability_test.mat', 'result', '-v7.3');
fprintf('\nSaved to matlab_realizability_test.mat\n');
