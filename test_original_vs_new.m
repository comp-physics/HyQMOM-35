% Test new implementation against original code
% Both should produce identical results

clc;
clear;
close all;

fprintf('========================================\n');
fprintf('Testing New vs Original Implementation\n');
fprintf('========================================\n\n');

% Test parameters
Np = 40;
tmax = 0.1;

%% Run new serial implementation
fprintf('Running NEW implementation (serial)...\n');
tic;
new_results = main(Np, tmax, false);
new_time = toc;
fprintf('  Completed in %.2f seconds\n', new_time);

%% Run original implementation (modified to Np=40)
fprintf('\nRunning ORIGINAL implementation...\n');
fprintf('  (Need to modify original/main_2Dcrossing_3DHyQMOM35.m to use Np=40)\n');
fprintf('  Skipping for now - will create modified version\n');

% Save new results for comparison
save('test_new_Np40_tmax01.mat', 'new_results', 'Np', 'tmax');
fprintf('\nSaved new results to test_new_Np40_tmax01.mat\n');

fprintf('\nNew implementation results:\n');
fprintf('  Grid: %dx%d\n', Np, Np);
fprintf('  Final time: %.4f\n', new_results.parameters.final_time);
fprintf('  Number of timesteps: %d\n', new_results.parameters.num_timesteps);
fprintf('  M(20,20,1) = %.10f\n', new_results.moments.M(20,20,1));
fprintf('  M(20,20,2) = %.10f\n', new_results.moments.M(20,20,2));

fprintf('\n========================================\n');
fprintf('To compare with original:\n');
fprintf('1. Modify original/main_2Dcrossing_3DHyQMOM35.m:\n');
fprintf('   - Set Np = 40 (line 30)\n');
fprintf('   - Set tmax = 0.1 (line 24)\n');
fprintf('2. Run original code\n');
fprintf('3. Compare M arrays\n');
fprintf('========================================\n');

