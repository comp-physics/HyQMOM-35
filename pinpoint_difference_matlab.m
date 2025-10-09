% Pinpoint exactly where M210 and M130 become nonzero in MATLAB
% Track cell (8,9) through realizability corrections

clear all;
close all;

setup_paths;

fprintf('===== PINPOINTING DIFFERENCE - MATLAB =====\n\n');

% Same setup as simulation
Np = 20;
rho = 0.7458726689777664;  % From cell (8,9) at step 1
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

fprintf('Testing InitializeM4_35 with actual simulation values:\n');
fprintf('  rho = %.15e\n', rho);
fprintf('  umean = %.15e\n', umean);
fprintf('  C200 = %.15e\n', C200);
fprintf('  C110 = %.15e (should be 0)\n\n', C110);

% Call InitializeM4_35
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);

fprintf('Output from InitializeM4_35:\n');
fprintf('  M210 = %.15e\n', M(8));
fprintf('  M130 = %.15e\n\n', M(14));

% Now apply realizability corrections
fprintf('=== APPLYING REALIZABILITY ===\n');

% First Flux_closure35
[~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(M, flag2D, Ma);
fprintf('After Flux_closure35_and_realizable_3D:\n');
fprintf('  M210 = %.15e\n', Mr(8));
fprintf('  M130 = %.15e\n\n', Mr(14));

% Eigenvalues X
[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
fprintf('After eigenvalues6_hyperbolic_3D (x-direction):\n');
fprintf('  M210 = %.15e\n', Mr(8));
fprintf('  M130 = %.15e\n\n', Mr(14));

% Eigenvalues Y  
[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
fprintf('After eigenvalues6_hyperbolic_3D (y-direction):\n');
fprintf('  M210 = %.15e\n', Mr(8));
fprintf('  M130 = %.15e\n\n', Mr(14));

% Second Flux_closure35
[~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
fprintf('After second Flux_closure35_and_realizable_3D:\n');
fprintf('  M210 = %.15e\n', Mr(8));
fprintf('  M130 = %.15e\n\n', Mr(14));

% Apply collision
dt = 0.02886751345948128;
Kn = 1.0;

fprintf('=== BEFORE CALLING COLLISION35 ===\n');
fprintf('  Mr(1) = %.15e\n', Mr(1));
fprintf('  Mr(2) = %.15e\n', Mr(2));
fprintf('  Mr(3) = %.15e\n', Mr(3));
fprintf('  Mr(6) = %.15e\n', Mr(6));
fprintf('  Mr(8) = %.15e\n', Mr(8));
fprintf('  Mr(10) = %.15e\n', Mr(10));
fprintf('  Mr(14) = %.15e\n', Mr(14));
fprintf('  Mr(16) = %.15e\n', Mr(16));
fprintf('  Mr(20) = %.15e\n\n', Mr(20));

Mr_final = collision35(Mr, dt, Kn);

fprintf('After collision35:\n');
fprintf('  M210 = %.15e\n', Mr_final(8));
fprintf('  M130 = %.15e\n\n', Mr_final(14));

fprintf('=== SUMMARY ===\n');
fprintf('Initial (from InitializeM4_35):  M210 = %.15e, M130 = %.15e\n', M(8), M(14));
fprintf('Final (after all corrections):   M210 = %.15e, M130 = %.15e\n', Mr_final(8), Mr_final(14));

fprintf('\n===== END MATLAB =====\n');
