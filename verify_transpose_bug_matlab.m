% Verify the transpose bug

clear all;
close all;

setup_paths;

fprintf('===== VERIFYING TRANSPOSE BUG - MATLAB =====\n\n');

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
[~, ~, ~, Mr1] = Flux_closure35_and_realizable_3D(M1, flag2D, Ma);

fprintf('M1 is %s\n', mat2str(size(M1)));
fprintf('Mr1 is %s\n\n', mat2str(size(Mr1)));

fprintf('=== TEST 1: collision35(Mr1) - row vector ===\n');
Mout_row = collision35(Mr1, dt, Kn);
fprintf('Result: M(8) = %.15e (WRONG: %.3f)\n\n', Mout_row(8), Mout_row(8));

fprintf('=== TEST 2: collision35(Mr1'') - transposed to column ===\n');
Mout_col = collision35(Mr1', dt, Kn);
fprintf('Result: M(8) = %.15e (should be 0)\n\n', Mout_col(8));

fprintf('=== VERIFICATION ===\n');
if abs(Mout_col(8)) < 0.01
    fprintf('✅ BUG CONFIRMED: Transpose fixes the issue!\n');
    fprintf('   Row vector gives: %.3f\n', Mout_row(8));
    fprintf('   Column vector gives: %.6f\n', Mout_col(8));
else
    fprintf('❌ Transpose doesn''t fix it\n');
end

fprintf('\n=== CHECKING WHAT collision35 SEES ===\n');
fprintf('When Mr1 is a row vector [1x35]:\n');
fprintf('  Mr1(1) = %.15e (rho)\n', Mr1(1));
fprintf('  Mr1(2) = %.15e (should be rho*u)\n', Mr1(2));
fprintf('  Mr1(3) = %.15e\n', Mr1(3));
fprintf('  Mr1(6) = %.15e (should be rho*v)\n', Mr1(6));
fprintf('  Mr1(8) = %.15e (M210)\n', Mr1(8));
fprintf('  Mr1(10) = %.15e\n', Mr1(10));
fprintf('  Mr1(16) = %.15e (should be rho*w)\n', Mr1(16));
fprintf('  Mr1(20) = %.15e\n\n', Mr1(20));

fprintf('Computing Theta from row vector:\n');
rho_from_row = Mr1(1);
umean_from_row = Mr1(2) / rho_from_row;
vmean_from_row = Mr1(6) / rho_from_row;
wmean_from_row = Mr1(16) / rho_from_row;
C200_from_row = Mr1(3)/rho_from_row - umean_from_row^2;
C020_from_row = Mr1(10)/rho_from_row - vmean_from_row^2;
C002_from_row = Mr1(20)/rho_from_row - wmean_from_row^2;
Theta_from_row = (C200_from_row + C020_from_row + C002_from_row) / 3;
fprintf('  Theta = %.15e\n', Theta_from_row);
fprintf('  Wrong output M(8) = %.15e\n', Mout_row(8));
fprintf('  Ratio = %.6f\n', Mout_row(8) / Theta_from_row);

fprintf('\n===== END =====\n');
