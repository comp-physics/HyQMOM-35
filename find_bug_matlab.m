% Systematically find the bug in MATLAB collision operator

clear all;
close all;

setup_paths;

fprintf('===== FINDING THE BUG - MATLAB =====\n\n');

% Same input as pinpoint_difference_matlab.m
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

fprintf('=== TEST 1: Direct call to collision35 ===\n');
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
fprintf('Before collision: M(8)=%.15e, M(14)=%.15e\n', M(8), M(14));
Mout1 = collision35(M, dt, Kn);
fprintf('After collision:  M(8)=%.15e, M(14)=%.15e\n\n', Mout1(8), Mout1(14));

fprintf('=== TEST 2: After one realizability call ===\n');
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
[~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(M, flag2D, Ma);
fprintf('After Flux_closure35: M(8)=%.15e, M(14)=%.15e\n', Mr(8), Mr(14));
Mout2 = collision35(Mr, dt, Kn);
fprintf('After collision:      M(8)=%.15e, M(14)=%.15e\n\n', Mout2(8), Mout2(14));

fprintf('=== TEST 3: After eigenvalues6_hyperbolic_3D (x) ===\n');
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
[~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(M, flag2D, Ma);
[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
fprintf('After eigenvalues6 (x): M(8)=%.15e, M(14)=%.15e\n', Mr(8), Mr(14));
Mout3 = collision35(Mr, dt, Kn);
fprintf('After collision:        M(8)=%.15e, M(14)=%.15e\n\n', Mout3(8), Mout3(14));

fprintf('=== TEST 4: After eigenvalues6_hyperbolic_3D (x and y) ===\n');
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
[~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(M, flag2D, Ma);
[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
fprintf('After eigenvalues6 (x,y): M(8)=%.15e, M(14)=%.15e\n', Mr(8), Mr(14));
Mout4 = collision35(Mr, dt, Kn);
fprintf('After collision:          M(8)=%.15e, M(14)=%.15e\n\n', Mout4(8), Mout4(14));

fprintf('=== TEST 5: Full sequence (as in pinpoint_difference) ===\n');
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
[~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(M, flag2D, Ma);
[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
[~, ~, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
[~, ~, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
fprintf('Before collision: M(8)=%.15e, M(14)=%.15e\n', Mr(8), Mr(14));
fprintf('  M(1)=%.15e, M(2)=%.15e, M(3)=%.15e\n', Mr(1), Mr(2), Mr(3));
fprintf('  M(6)=%.15e, M(10)=%.15e, M(16)=%.15e, M(20)=%.15e\n', Mr(6), Mr(10), Mr(16), Mr(20));

% Create a copy to ensure no aliasing
Mr_copy = Mr;
fprintf('\nCalling collision35 with copy...\n');
Mout5 = collision35(Mr_copy, dt, Kn);
fprintf('After collision:  M(8)=%.15e, M(14)=%.15e\n\n', Mout5(8), Mout5(14));

% Check if Mr was modified
if any(Mr ~= Mr_copy)
    fprintf('WARNING: Mr was modified by collision35!\n');
end

fprintf('=== CHECKING WHICH TEST GIVES WRONG ANSWER ===\n');
if abs(Mout1(8)) > 0.1
    fprintf('Test 1 (direct): WRONG (%.3f)\n', Mout1(8));
else
    fprintf('Test 1 (direct): correct (0.0)\n');
end

if abs(Mout2(8)) > 0.1
    fprintf('Test 2 (after Flux_closure35): WRONG (%.3f)\n', Mout2(8));
else
    fprintf('Test 2 (after Flux_closure35): correct (0.0)\n');
end

if abs(Mout3(8)) > 0.1
    fprintf('Test 3 (after eigenvalues x): WRONG (%.3f)\n', Mout3(8));
else
    fprintf('Test 3 (after eigenvalues x): correct (0.0)\n');
end

if abs(Mout4(8)) > 0.1
    fprintf('Test 4 (after eigenvalues x,y): WRONG (%.3f)\n', Mout4(8));
else
    fprintf('Test 4 (after eigenvalues x,y): correct (0.0)\n');
end

if abs(Mout5(8)) > 0.1
    fprintf('Test 5 (full sequence): WRONG (%.3f)\n', Mout5(8));
else
    fprintf('Test 5 (full sequence): correct (0.0)\n');
end

fprintf('\n===== END =====\n');
