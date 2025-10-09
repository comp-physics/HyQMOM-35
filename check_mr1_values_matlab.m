% Check actual values in Mr1

clear all;
close all;

setup_paths;

fprintf('===== CHECKING Mr1 VALUES - MATLAB =====\n\n');

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

M1 = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
fprintf('M1 from InitializeM4_35:\n');
fprintf('  M1(1)  = %.15e\n', M1(1));
fprintf('  M1(2)  = %.15e\n', M1(2));
fprintf('  M1(3)  = %.15e\n', M1(3));
fprintf('  M1(6)  = %.15e\n', M1(6));
fprintf('  M1(8)  = %.15e\n', M1(8));
fprintf('  M1(10) = %.15e\n', M1(10));
fprintf('  M1(14) = %.15e\n', M1(14));
fprintf('  M1(16) = %.15e\n', M1(16));
fprintf('  M1(20) = %.15e\n\n', M1(20));

[~, ~, ~, Mr1] = Flux_closure35_and_realizable_3D(M1, flag2D, Ma);
fprintf('Mr1 from Flux_closure35_and_realizable_3D:\n');
fprintf('  Mr1(1)  = %.15e\n', Mr1(1));
fprintf('  Mr1(2)  = %.15e\n', Mr1(2));
fprintf('  Mr1(3)  = %.15e\n', Mr1(3));
fprintf('  Mr1(6)  = %.15e\n', Mr1(6));
fprintf('  Mr1(8)  = %.15e\n', Mr1(8));
fprintf('  Mr1(10) = %.15e\n', Mr1(10));
fprintf('  Mr1(14) = %.15e\n', Mr1(14));
fprintf('  Mr1(16) = %.15e\n', Mr1(16));
fprintf('  Mr1(20) = %.15e\n\n', Mr1(20));

fprintf('Comparison:\n');
fprintf('  M1(8) - Mr1(8) = %.15e\n', M1(8) - Mr1(8));
fprintf('  M1(14) - Mr1(14) = %.15e\n\n', M1(14) - Mr1(14));

fprintf('Are they identical?\n');
if isequal(M1, Mr1)
    fprintf('  YES - M1 and Mr1 are identical\n');
else
    fprintf('  NO - M1 and Mr1 differ\n');
    fprintf('  Differences at indices:\n');
    diffs = find(abs(M1 - Mr1) > 1e-14);
    for idx = diffs'
        fprintf('    Index %d: M1=%.15e, Mr1=%.15e, diff=%.15e\n', ...
                idx, M1(idx), Mr1(idx), M1(idx)-Mr1(idx));
    end
end

fprintf('\n=== NOW CALL collision35 ON BOTH ===\n');
Mout1 = collision35(M1, dt, Kn);
Mout_r1 = collision35(Mr1, dt, Kn);

fprintf('collision35(M1):\n');
fprintf('  M(8) = %.15e\n', Mout1(8));
fprintf('collision35(Mr1):\n');
fprintf('  M(8) = %.15e\n\n', Mout_r1(8));

fprintf('===== END =====\n');
