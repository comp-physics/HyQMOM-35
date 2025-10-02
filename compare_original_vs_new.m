% Comprehensive comparison: Original vs New Implementation
% Both with Np=40, tmax=0.1

clc;
clear;
close all;

fprintf('============================================================\n');
fprintf('   Comparing Original vs New Implementation\n');
fprintf('   Both using Np=40, tmax=0.1\n');
fprintf('============================================================\n\n');

%% Run NEW implementation
fprintf('1. Running NEW implementation (serial)...\n');
tic;
new_results = main(40, 0.1, false);
new_time = toc;
fprintf('   ✓ Completed in %.2f seconds\n', new_time);
if isfield(new_results, 'parameters')
    if isfield(new_results.parameters, 'num_timesteps')
        fprintf('   - Timesteps: %d\n', new_results.parameters.num_timesteps);
    end
    if isfield(new_results.parameters, 'final_time')
        fprintf('   - Final time: %.6f\n', new_results.parameters.final_time);
    end
end

%% Run ORIGINAL implementation  
fprintf('\n2. Running ORIGINAL implementation...\n');
% Need to run in original directory with proper path
current_dir = pwd;
cd('original');
% Add original directory to path temporarily
addpath(pwd);
tic;
run('../test_original_Np40.m');
cd(current_dir);
original_time = toc;
fprintf('   ✓ Completed in %.2f seconds\n', original_time);

%% Load results
fprintf('\n3. Loading and comparing results...\n');
load('test_original_Np40_results.mat', 'original_results');

%% Compare
fprintf('\n============================================================\n');
fprintf('   COMPARISON RESULTS\n');
fprintf('============================================================\n\n');

% Compare M arrays
M_new = new_results.moments.M;
M_orig = original_results.M;

diff_M = abs(M_new - M_orig);
max_diff_M = max(diff_M(:));
rel_diff_M = max_diff_M / max(abs(M_orig(:)));

fprintf('Moment Array (M):\n');
fprintf('  Max absolute difference: %.15e\n', max_diff_M);
fprintf('  Max relative difference: %.15e\n', rel_diff_M);
fprintf('  New    M(20,20,1): %.15e\n', M_new(20,20,1));
fprintf('  Original M(20,20,1): %.15e\n', M_orig(20,20,1));
fprintf('  Difference:          %.15e\n', abs(M_new(20,20,1) - M_orig(20,20,1)));

% Find location of max difference
[~, idx] = max(diff_M(:));
[i, j, k] = ind2sub(size(M_new), idx);
fprintf('\n  Max difference at M(%d,%d,%d):\n', i, j, k);
fprintf('    New:      %.15e\n', M_new(i,j,k));
fprintf('    Original: %.15e\n', M_orig(i,j,k));
fprintf('    Diff:     %.15e\n', diff_M(i,j,k));

% Compare C arrays
if isfield(new_results.moments, 'C')
    C_new = new_results.moments.C;
    C_orig = original_results.C;
    
    diff_C = abs(C_new - C_orig);
    max_diff_C = max(diff_C(:));
    
    fprintf('\nCumulant Array (C):\n');
    fprintf('  Max absolute difference: %.15e\n', max_diff_C);
end

% Compare S arrays
if isfield(new_results.moments, 'S')
    S_new = new_results.moments.S;
    S_orig = original_results.S;
    
    diff_S = abs(S_new - S_orig);
    max_diff_S = max(diff_S(:));
    
    fprintf('\nSecond-order moments (S):\n');
    fprintf('  Max absolute difference: %.15e\n', max_diff_S);
end

% Overall assessment
fprintf('\n============================================================\n');
fprintf('   ASSESSMENT\n');
fprintf('============================================================\n\n');

if max_diff_M < 1e-14
    fprintf('✅ PERFECT MATCH - Bitwise identical (within machine precision)\n');
elseif max_diff_M < 1e-10
    fprintf('✅ EXCELLENT - Differences at floating-point roundoff level\n');
elseif max_diff_M < 1e-6
    fprintf('✅ VERY GOOD - Differences at expected numerical precision\n');
elseif max_diff_M < 1e-3
    fprintf('⚠️  ACCEPTABLE - Small differences, may need investigation\n');
else
    fprintf('❌ SIGNIFICANT DIFFERENCES - Implementation may differ\n');
end

fprintf('\nMax difference: %.6e (relative: %.6e)\n', max_diff_M, rel_diff_M);
fprintf('\nTiming:\n');
fprintf('  Original: %.2f seconds\n', original_time);
fprintf('  New:      %.2f seconds\n', new_time);
fprintf('  Speedup:  %.2fx\n', original_time/new_time);

fprintf('\n============================================================\n');

