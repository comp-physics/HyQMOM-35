% Simplified debug: just run simulation_runner and print values at each step
% This matches exactly what MATLAB does

clear; clc;

% Add paths
addpath('src');
addpath('src/autogen');

% Parameters
params = struct();
params.Np = 20;
params.tmax = 0.1;
params.Kn = 1.0;
params.Ma = 0.0;
params.flag2D = 0;
params.CFL = 0.5;
params.dx = 1.0 / params.Np;
params.dy = 1.0 / params.Np;
params.Nmom = 35;
params.nnmax = 1000;
params.dtmax = params.CFL * min(params.dx, params.dy);
params.rhol = 1.0;
params.rhor = 0.01;
params.T = 1.0;
params.r110 = 0.0;
params.r101 = 0.0;
params.r011 = 0.0;
params.symmetry_check_interval = 10;
params.enable_memory_tracking = false;

% Modify simulation_runner to add debug output
% For now, just run it normally and check if it completes

fprintf('Running simulation with Np=%d, tmax=%.2f\n', params.Np, params.tmax);
fprintf('Expected: 13 steps to reach tmax\n\n');

% Run simulation
try
    results = simulation_runner(params);
    fprintf('\nSimulation completed successfully!\n');
    fprintf('Final time: %.15e\n', results.t);
    fprintf('Number of steps: %d\n', results.nn);
    fprintf('Final M(1,1,1:5) = [%.15e, %.15e, %.15e, %.15e, %.15e]\n', ...
        results.M(1,1,1), results.M(1,1,2), results.M(1,1,3), results.M(1,1,4), results.M(1,1,5));
catch ME
    fprintf('\nSimulation FAILED with error:\n');
    fprintf('%s\n', ME.message);
    fprintf('Error occurred in: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end
