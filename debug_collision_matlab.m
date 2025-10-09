% Debug collision operator step by step

clear all;
close all;

setup_paths;

fprintf('===== DEBUG COLLISION - MATLAB =====\n\n');

% Input moments
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

fprintf('=== INPUT TO COLLISION ===\n');
M = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
fprintf('M(1) = %.15e (rho)\n', M(1));
fprintf('M(2) = %.15e (rho*u)\n', M(2));
fprintf('M(3) = %.15e\n', M(3));
fprintf('M(6) = %.15e (rho*v)\n', M(6));
fprintf('M(8) = %.15e (M210)\n', M(8));
fprintf('M(10) = %.15e\n', M(10));
fprintf('M(14) = %.15e (M130)\n', M(14));
fprintf('M(16) = %.15e (rho*w)\n', M(16));
fprintf('M(20) = %.15e\n\n', M(20));

% Collision parameters
dt = 0.02886751345948128;
Kn = 1.0;

fprintf('=== INSIDE COLLISION OPERATOR ===\n');

% Extract quantities (same as collision35)
rho_extracted = M(1);
umean_extracted = M(2) / rho_extracted;
vmean_extracted = M(6) / rho_extracted;
wmean_extracted = M(16) / rho_extracted;

fprintf('Extracted values:\n');
fprintf('  rho = %.15e\n', rho_extracted);
fprintf('  umean = %.15e\n', umean_extracted);
fprintf('  vmean = %.15e\n', vmean_extracted);
fprintf('  wmean = %.15e\n\n', wmean_extracted);

% Compute temperature
C200_extracted = M(3)/rho_extracted - umean_extracted^2;
C020_extracted = M(10)/rho_extracted - vmean_extracted^2;
C002_extracted = M(20)/rho_extracted - wmean_extracted^2;
Theta = (C200_extracted + C020_extracted + C002_extracted) / 3;

fprintf('Covariance components:\n');
fprintf('  C200 = %.15e\n', C200_extracted);
fprintf('  C020 = %.15e\n', C020_extracted);
fprintf('  C002 = %.15e\n', C002_extracted);
fprintf('  Theta = %.15e\n\n', Theta);

% Collision time scale
tc = Kn / (rho_extracted * sqrt(Theta) * 2);
fprintf('Collision time scale:\n');
fprintf('  tc = %.15e\n', tc);
fprintf('  exp(-dt/tc) = %.15e\n\n', exp(-dt/tc));

% Maxwellian equilibrium
fprintf('=== MAXWELLIAN EQUILIBRIUM ===\n');
MG = InitializeM4_35(rho_extracted, umean_extracted, vmean_extracted, wmean_extracted, Theta, 0, 0, Theta, 0, Theta);
fprintf('MG(1) = %.15e\n', MG(1));
fprintf('MG(2) = %.15e\n', MG(2));
fprintf('MG(3) = %.15e\n', MG(3));
fprintf('MG(8) = %.15e (M210)\n', MG(8));
fprintf('MG(14) = %.15e (M130)\n\n', MG(14));

% BGK relaxation
fprintf('=== BGK RELAXATION ===\n');
Mout = MG - exp(-dt/tc) * (MG - M);
fprintf('Mout(1) = %.15e\n', Mout(1));
fprintf('Mout(2) = %.15e\n', Mout(2));
fprintf('Mout(3) = %.15e\n', Mout(3));
fprintf('Mout(8) = %.15e (M210)\n', Mout(8));
fprintf('Mout(14) = %.15e (M130)\n\n', Mout(14));

% Check specific moment
fprintf('=== DETAILED M210 CALCULATION ===\n');
fprintf('  MG(8)  = %.15e\n', MG(8));
fprintf('  M(8)   = %.15e\n', M(8));
fprintf('  MG(8) - M(8) = %.15e\n', MG(8) - M(8));
fprintf('  exp(-dt/tc) * (MG(8) - M(8)) = %.15e\n', exp(-dt/tc) * (MG(8) - M(8)));
fprintf('  Mout(8) = MG(8) - exp(-dt/tc) * (MG(8) - M(8)) = %.15e\n', Mout(8));

fprintf('\n===== END MATLAB DEBUG =====\n');
