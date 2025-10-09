% Check sizes of M1 and Mr1

clear all;
close all;

setup_paths;

fprintf('===== CHECKING SIZES - MATLAB =====\n\n');

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

M1 = InitializeM4_35(rho, umean, vmean, wmean, C200, C110, C101, C020, C011, C002);
fprintf('M1:\n');
fprintf('  size(M1) = %s\n', mat2str(size(M1)));
fprintf('  length(M1) = %d\n', length(M1));
fprintf('  class(M1) = %s\n', class(M1));
fprintf('  M1(end) = %.15e\n\n', M1(end));

[~, ~, ~, Mr1] = Flux_closure35_and_realizable_3D(M1, flag2D, Ma);
fprintf('Mr1:\n');
fprintf('  size(Mr1) = %s\n', mat2str(size(Mr1)));
fprintf('  length(Mr1) = %d\n', length(Mr1));
fprintf('  class(Mr1) = %s\n', class(Mr1));
fprintf('  Mr1(end) = %.15e\n\n', Mr1(end));

if length(Mr1) > 35
    fprintf('❌ BUG FOUND: Mr1 has %d elements instead of 35!\n', length(Mr1));
    fprintf('Extra elements:\n');
    for i = 36:length(Mr1)
        fprintf('  Mr1(%d) = %.15e\n', i, Mr1(i));
    end
end

fprintf('\n===== END =====\n');
