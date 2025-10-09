% Investigate why MATLAB produces ~1e-9 where Julia produces 0
% Focus on InitializeM4_35 and S4toC4_3D_r

clear all;
close all;

setup_paths;

fprintf('===== INVESTIGATING ZERO vs 1e-9 - MATLAB =====\n\n');

% Use exact values from cell (8,9) at step 1
rho = 0.7458726689777664;
umean = -0.13731136169066;
vmean = 0.0;
wmean = 0.0;
C200 = 0.6829490116088;
C020 = 0.7458726689777664;
C002 = 0.7458726689777664;
C110 = 0.0;  % Exactly zero - no correlation
C101 = 0.0;
C011 = 0.0;

fprintf('=== INPUT PARAMETERS ===\n');
fprintf('Covariance matrix inputs:\n');
fprintf('  C200 = %.15e\n', C200);
fprintf('  C110 = %.15e (exactly zero)\n', C110);
fprintf('  C101 = %.15e (exactly zero)\n', C101);
fprintf('  C020 = %.15e\n', C020);
fprintf('  C011 = %.15e (exactly zero)\n', C011);
fprintf('  C002 = %.15e\n\n', C002);

% Step 1: Call S4toC4_3D_r directly
fprintf('=== STEP 1: S4toC4_3D_r ===\n');

% Standardized moments for Gaussian
S300=0; S210=0; S201=0; S120=0; S111=0; S102=0; S030=0; S021=0; S012=0; S003=0;
S400=3; S310=0; S301=0; S220=1; S211=0; S202=1;
S130=0; S121=0; S112=0; S103=0; S040=3; S031=0; S022=1; S013=0; S004=3;

fprintf('Calling S4toC4_3D_r...\n');
C4_output = S4toC4_3D_r(C200,C110,C101,C020,C011,C002,...
                        S300,S210,S201,S120,S111,S102,S030,S021,S012,S003,...
                        S400,S310,S301,S220,S211,S202,S130,S121,S112,S103,S040,S031,S022,S013,S004);

fprintf('C4_output size: %s\n', mat2str(size(C4_output)));
fprintf('C4_output(8) = %.15e (this becomes M210)\n', C4_output(8));
fprintf('C4_output(14) = %.15e (this becomes M130)\n\n', C4_output(14));

% Step 2: Look at matrix square root in S4toC4_3D_r
fprintf('=== STEP 2: Matrix Square Root Analysis ===\n');

% Build the covariance matrix
C2 = [C200 C110 C101; C110 C020 C011; C101 C011 C002];
fprintf('Covariance matrix C2:\n');
disp(C2);

fprintf('\nComputing sqrtm(C2)...\n');
A = sqrtm(C2);
fprintf('Matrix square root A:\n');
disp(A);

fprintf('\nChecking A*A - C2:\n');
diff = A*A - C2;
disp(diff);
fprintf('Max difference: %.15e\n\n', max(abs(diff(:))));

% Step 3: Check eigenvalues and eigenvectors
fprintf('=== STEP 3: Eigenvalue Analysis ===\n');
[V, D] = eig(C2);
fprintf('Eigenvalues of C2:\n');
disp(diag(D));
fprintf('\nEigenvectors of C2:\n');
disp(V);

% Step 4: Trace through InitializeM4_35
fprintf('\n=== STEP 4: InitializeM4_35 Output ===\n');
M_output = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);

fprintf('M_output size: %s\n', mat2str(size(M_output)));
fprintf('M_output(8)  = %.15e (M210)\n', M_output(8));
fprintf('M_output(14) = %.15e (M130)\n\n', M_output(14));

% Step 5: Check if it's a precision/tolerance issue
fprintf('=== STEP 5: Precision Check ===\n');
fprintf('eps = %.15e\n', eps);
fprintf('realmin = %.15e\n', realmin);
fprintf('Is M(8) exactly zero? %s\n', mat2str(M_output(8) == 0));
fprintf('Is M(8) < eps? %s\n', mat2str(abs(M_output(8)) < eps));
fprintf('Is M(8) < 1e-12? %s\n\n', mat2str(abs(M_output(8)) < 1e-12));

% Step 6: Manual calculation to understand the computation
fprintf('=== STEP 6: Manual Calculation ===\n');
fprintf('From central to raw moments:\n');

% C210 is the central moment we want
fprintf('Looking for C210 (central moment)...\n');

% The formula involves sqrt(C200) * sqrt(C020) terms
fprintf('sqrt(C200) = %.15e\n', sqrt(C200));
fprintf('sqrt(C020) = %.15e\n', sqrt(C020));
fprintf('sqrt(C002) = %.15e\n', sqrt(C002));

% Extract from C4_output
[C000, C100, C200_out, C300, C400, C010, C110_out, C210, C310, ...
 C020_out, C120, C220, C030, C130, C040] = moment_conversion_utils('M4_to_vars', C4_output);

fprintf('\nExtracted central moments:\n');
fprintf('  C210 = %.15e\n', C210);
fprintf('  C120 = %.15e\n', C120);
fprintf('  C130 = %.15e\n', C130);

fprintf('\nComputing raw moment M210 from central C210:\n');
fprintf('  M210 = C210 + umean*C020 + 2*vmean*C110 + vmean^2*C100\n');
fprintf('       + umean^2*C010 + 2*umean*vmean*C010 + vmean^2*umean*M000\n');

M210_manual = C210 + umean*C020_out + 2*vmean*C110_out;
fprintf('  Simplified (vmean=0): M210 = C210 + umean*C020\n');
fprintf('  M210_manual = %.15e + %.15e*%.15e\n', C210, umean, C020_out);
fprintf('  M210_manual = %.15e\n\n', M210_manual);

fprintf('=== SUMMARY ===\n');
fprintf('MATLAB produces M210 = %.15e\n', M_output(8));
fprintf('This is %s zero\n', iif(M_output(8)==0, 'EXACTLY', 'NOT'));
fprintf('Magnitude: %.3e\n', abs(M_output(8)));

fprintf('\n===== END MATLAB INVESTIGATION =====\n');

function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
