% Isolate the exact bug - step by step

clear all;
close all;
clearvars -global;  % Clear any global variables

setup_paths;

fprintf('===== ISOLATING BUG - MATLAB =====\n\n');

% Test parameters
rho = 0.7458726689777664;
umean = -0.13731136169066;
vmean = 0.0;
wmean = 0.0;
C200 = 0.6829490116088;
C020 = 0.7458726689777664;
C002 = 0.7458726689777664;
C110 = 0.0;
C101 = 0.0;
C011 = 0.0;
Ma = 0.0;
flag2D = 0;
dt = 0.02886751345948128;
Kn = 1.0;

fprintf('=== BEFORE ANY FUNCTION CALLS ===\n');
M_test = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
fprintf('Test collision35 (should be 0):\n');
Mout_before = collision35(M_test, dt, Kn);
fprintf('  Result: M(8)=%.15e\n\n', Mout_before(8));

fprintf('=== CALL Flux_closure35_and_realizable_3D ===\n');
M1 = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
[~, ~, ~, Mr1] = Flux_closure35_and_realizable_3D(M1, flag2D, Ma);
fprintf('Returned from Flux_closure35_and_realizable_3D\n');
fprintf('  Mr1(8) = %.15e (should be 0)\n\n', Mr1(8));

fprintf('=== NOW TEST collision35 AGAIN ===\n');
M_test2 = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
fprintf('Test collision35 (should still be 0):\n');
Mout_after = collision35(M_test2, dt, Kn);
fprintf('  Result: M(8)=%.15e\n', Mout_after(8));

if abs(Mout_after(8)) > 0.1
    fprintf('\n❌ BUG CONFIRMED: collision35 corrupted after Flux_closure35 call!\n');
    fprintf('   Value changed from %.3f to %.3f\n', Mout_before(8), Mout_after(8));
else
    fprintf('\n✅ No bug: collision35 still works correctly\n');
end

fprintf('\n=== CHECK: Call collision35 on Mr1 directly ===\n');
Mout_mr1 = collision35(Mr1, dt, Kn);
fprintf('  collision35(Mr1): M(8)=%.15e\n', Mout_mr1(8));

fprintf('\n=== HYPOTHESIS: InitializeM4_35 was corrupted ===\n');
fprintf('Let me call InitializeM4_35 again and check output...\n');
M_new = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
fprintf('  M_new(8) = %.15e (should be 0)\n', M_new(8));

% Check if Theta is being stored somewhere
fprintf('\n=== DEBUG: Check Theta value ===\n');
Theta_computed = (C200 + C020 + C002) / 3;
fprintf('  Computed Theta = %.15e\n', Theta_computed);
fprintf('  Wrong M(8) value = %.15e\n', Mout_after(8));
fprintf('  Ratio = %.6f\n', Mout_after(8) / Theta_computed);

fprintf('\n===== END =====\n');
