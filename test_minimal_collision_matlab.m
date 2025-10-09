% Minimal test case: Collision operator on isotropic Maxwellian
% This should produce exactly zero for cross-moments M210 and M130
% but MATLAB produces ~1e-9 residuals

clear all;
close all;

% Add paths
setup_paths;

fprintf('===== MINIMAL COLLISION TEST - MATLAB =====\n\n');

% Test Case 1: Isotropic Maxwellian (all correlations = 0)
fprintf('Test Case 1: Isotropic Maxwellian\n');
fprintf('%s\n', repmat('-', 1, 50));

% Initial conditions (isotropic, at rest)
rho = 1.0;
umean = 0.0;
vmean = 0.0;
wmean = 0.0;
Theta = 1.0;  % Temperature

% Create initial Maxwellian with isotropic covariance
M_initial = InitializeM4_35(rho, umean, vmean, wmean, Theta, 0.0, 0.0, Theta, 0.0, Theta);

fprintf('Initial moments (should be isotropic Maxwellian):\n');
fprintf('  M000 = %.15e\n', M_initial(1));
fprintf('  M100 = %.15e\n', M_initial(2));
fprintf('  M200 = %.15e\n', M_initial(3));
fprintf('  M210 = %.15e  <-- Cross-moment (should be ~0)\n', M_initial(8));
fprintf('  M130 = %.15e  <-- Cross-moment (should be ~0)\n', M_initial(14));
fprintf('\n');

% Apply collision operator (should preserve Maxwellian)
dt = 0.01;
Kn = 1.0;
M_after = collision35(M_initial, dt, Kn);

fprintf('After collision operator:\n');
fprintf('  M000 = %.15e\n', M_after(1));
fprintf('  M100 = %.15e\n', M_after(2));
fprintf('  M200 = %.15e\n', M_after(3));
fprintf('  M210 = %.15e  <-- Cross-moment\n', M_after(8));
fprintf('  M130 = %.15e  <-- Cross-moment\n', M_after(14));
fprintf('\n');

fprintf('Changes:\n');
fprintf('  ΔM210 = %.15e\n', M_after(8) - M_initial(8));
fprintf('  ΔM130 = %.15e\n', M_after(14) - M_initial(14));
fprintf('\n');

% Test Case 2: Apply S4toC4_3D_r directly
fprintf('\nTest Case 2: Direct S4toC4_3D_r call\n');
fprintf('%s\n', repmat('-', 1, 50));

% Standardized moments for Gaussian
S300=0; S210=0; S201=0; S120=0; S111=0; S102=0; S030=0; S021=0; S012=0; S003=0;
S400=3; S310=0; S301=0; S220=1; S211=0; S202=1;
S130=0; S121=0; S112=0; S103=0; S040=3; S031=0; S022=1; S013=0; S004=3;

% Covariance (isotropic)
C200 = Theta;
C020 = Theta;
C002 = Theta;
C110 = 0.0;
C101 = 0.0;
C011 = 0.0;

fprintf('Input covariance matrix:\n');
fprintf('  C200 = %.15e, C110 = %.15e, C101 = %.15e\n', C200, C110, C101);
fprintf('  C020 = %.15e, C011 = %.15e\n', C020, C011);
fprintf('  C002 = %.15e\n', C002);
fprintf('\n');

% Call S4toC4_3D_r
C4 = S4toC4_3D_r(C200, C110, C101, C020, C011, C002, ...
                 S300, S210, S201, S120, S111, S102, S030, S021, S012, S003, ...
                 S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, S040, S031, S022, S013, S004);

% Extract specific central moments
[~, ~, ~, ~, ~, ~, C110_out, C210, ~, ~, C120, ~, C030, C130, ~] = ...
    moment_conversion_utils('M4_to_vars', C4);

fprintf('Output central moments:\n');
fprintf('  C110 = %.15e\n', C110_out);
fprintf('  C210 = %.15e  <-- Should be 0\n', C210);
fprintf('  C120 = %.15e  <-- Should be 0\n', C120);
fprintf('  C130 = %.15e  <-- Should be 0\n', C130);
fprintf('\n');

% Test Case 3: Matrix square root
fprintf('\nTest Case 3: Matrix square root behavior\n');
fprintf('%s\n', repmat('-', 1, 50));

C2 = [C200 C110 C101; C110 C020 C011; C101 C011 C002];
fprintf('Covariance matrix C2:\n');
disp(C2);

A = sqrtm(C2);
fprintf('Matrix square root A = sqrtm(C2):\n');
disp(A);

fprintf('Check: A*A (should equal C2):\n');
disp(A * A);

fprintf('Difference: A*A - C2:\n');
disp(A * A - C2);

fprintf('Max absolute error: %.15e\n', max(abs(A*A - C2), [], 'all'));
fprintf('\n');

% Test Case 4: Small perturbation
fprintf('\nTest Case 4: Effect of tiny perturbation\n');
fprintf('%s\n', repmat('-', 1, 50));

% Add a tiny cross-correlation
epsilon = 1e-9;
M_perturbed = InitializeM4_35(rho, umean, vmean, wmean, Theta, epsilon, 0.0, Theta, 0.0, Theta);

fprintf('With tiny perturbation (r110 = %.2e):\n', epsilon);
fprintf('  M210 = %.15e\n', M_perturbed(8));
fprintf('  M130 = %.15e\n', M_perturbed(14));
fprintf('\n');

% Apply collision
M_perturbed_after = collision35(M_perturbed, dt, Kn);

fprintf('After collision:\n');
fprintf('  M210 = %.15e\n', M_perturbed_after(8));
fprintf('  M130 = %.15e\n', M_perturbed_after(14));
fprintf('\n');

fprintf('===== END MATLAB TEST =====\n');
