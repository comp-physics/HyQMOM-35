% Quick comparison of first two time steps

clear all;
close all;

setup_paths;

fprintf('===== MATLAB: Two Time Steps =====\n\n');

params = struct('Np', 20, 'Kn', 1.0, 'Ma', 0.0, 'tmax', 0.05, 'flag2D', 0, 'CFL', 0.5);

% Run for ~2 steps (adjust tmax if needed)
tic;
result = simulation_runner(params);
elapsed = toc;

fprintf('Completed in %.2f seconds\n', elapsed);
fprintf('Time steps: %d\n', result.time_steps);
fprintf('Final time: %.6f\n', result.final_time);

% Save results
M_matlab = result.M;
save('matlab_two_steps.mat', 'M_matlab', 'result');

fprintf('\nSample moments at cell (8,9):\n');
fprintf('  M(8,9,1) = %.15e\n', M_matlab(8,9,1));
fprintf('  M(8,9,3) = %.15e\n', M_matlab(8,9,3));
fprintf('  M(8,9,8) = %.15e\n', M_matlab(8,9,8));
fprintf('  M(8,9,14) = %.15e\n', M_matlab(8,9,14));

fprintf('\n✅ MATLAB results saved\n');
