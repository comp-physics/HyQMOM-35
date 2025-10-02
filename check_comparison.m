% Just check if results match - assuming both have already run
clc;

fprintf('Loading saved results...\n');

% Check if files exist
if ~exist('test_original_Np40_results.mat', 'file')
    error('Original results not found. Run test_original_Np40.m first');
end

if ~exist('test_new_Np40_tmax01.mat', 'file')
    error('New results not found. Run test_original_vs_new.m first');
end

% Load results
load('test_original_Np40_results.mat');
load('test_new_Np40_tmax01.mat');

fprintf('\n=== COMPARISON: Original vs New (Np=40, tmax=0.1) ===\n\n');

% Compare M arrays
M_new = new_results.moments.M;
M_orig = original_results.M;

diff_M = abs(M_new - M_orig);
max_diff = max(diff_M(:));
rel_diff = max_diff / max(abs(M_orig(:)));

fprintf('Moment Array (M):\n');
fprintf('  Size: %dx%dx%d\n', size(M_new));
fprintf('  Max absolute difference: %.15e\n', max_diff);
fprintf('  Max relative difference: %.15e\n', rel_diff);

fprintf('\nSample values at (20,20,1):\n');
fprintf('  New:      %.15e\n', M_new(20,20,1));
fprintf('  Original: %.15e\n', M_orig(20,20,1));
fprintf('  Diff:     %.15e\n', abs(M_new(20,20,1) - M_orig(20,20,1)));

% Find max difference location
[~, idx] = max(diff_M(:));
[i, j, k] = ind2sub(size(M_new), idx);
fprintf('\nMax difference at M(%d,%d,%d):\n', i, j, k);
fprintf('  New:      %.15e\n', M_new(i,j,k));
fprintf('  Original: %.15e\n', M_orig(i,j,k));
fprintf('  Diff:     %.15e\n', diff_M(i,j,k));

fprintf('\n=== ASSESSMENT ===\n');
if max_diff < 1e-14
    fprintf('✅ PERFECT - Bitwise identical\n');
elseif max_diff < 1e-10
    fprintf('✅ EXCELLENT - Roundoff level\n');
elseif max_diff < 1e-6
    fprintf('✅ VERY GOOD - Numerical precision\n');
else
    fprintf('⚠️  Differences present: %.3e\n', max_diff);
end

fprintf('\n');

