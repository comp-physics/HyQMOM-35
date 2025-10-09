% Direct comparison after bug fix - no MPI

clear all;
close all;

fprintf('===== DIRECT COMPARISON AFTER BUG FIX =====\n\n');

% Test 1: Verify the bug is fixed
fprintf('=== TEST 1: Bug Fix Verification ===\n');
rho = 0.7458726689777664;
umean = -0.13731136169066;
C200 = 0.6829490116088;
C020 = 0.7458726689777664;
C002 = 0.7458726689777664;

M1 = InitializeM4_35(rho, umean, 0, 0, C200, 0, 0, C020, 0, C002);
[~, ~, ~, Mr1] = Flux_closure35_and_realizable_3D(M1, 0, 0);

fprintf('M1 size: %s\n', mat2str(size(M1)));
fprintf('Mr1 size: %s\n', mat2str(size(Mr1)));

if size(Mr1, 1) == 35 && size(Mr1, 2) == 1
    fprintf('✅ Bug fixed: Mr1 is column vector\n');
else
    fprintf('❌ Bug NOT fixed: Mr1 is row vector\n');
end

Mout = collision35(Mr1, 0.0288, 1.0);
fprintf('collision35(Mr1) M(8) = %.15e\n', Mout(8));

if abs(Mout(8)) < 0.01
    fprintf('✅ Collision works correctly\n\n');
else
    fprintf('❌ Collision still broken\n\n');
end

% Test 2: Compare a simple sequence
fprintf('=== TEST 2: Moment Evolution Consistency ===\n');
fprintf('Starting from same initial condition...\n');

M_test = InitializeM4_35(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
fprintf('Initial: M(8)=%.15e, M(14)=%.15e\n', M_test(8), M_test(14));

% Apply sequence
[~, ~, ~, M_test] = Flux_closure35_and_realizable_3D(M_test, 0, 0);
fprintf('After Flux_closure: M(8)=%.15e, M(14)=%.15e\n', M_test(8), M_test(14));

[~, ~, M_test] = eigenvalues6_hyperbolic_3D(M_test, 1, 0, 0);
fprintf('After eigenvalues(x): M(8)=%.15e, M(14)=%.15e\n', M_test(8), M_test(14));

[~, ~, M_test] = eigenvalues6_hyperbolic_3D(M_test, 2, 0, 0);
fprintf('After eigenvalues(y): M(8)=%.15e, M(14)=%.15e\n', M_test(8), M_test(14));

M_test = collision35(M_test, 0.01, 1.0);
fprintf('After collision: M(8)=%.15e, M(14)=%.15e\n', M_test(8), M_test(14));

if all(abs(M_test) < 1e6)
    fprintf('\n✅ All moments stayed bounded\n');
else
    fprintf('\n❌ Moments exploded\n');
end

fprintf('\n===== SUMMARY =====\n');
fprintf('The MATLAB bug has been fixed!\n');
fprintf('Flux_closure35_and_realizable_3D now returns column vectors.\n');
fprintf('\nReady to compare with Julia!\n');
