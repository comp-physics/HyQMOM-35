function create_test_goldenfiles()
%CREATE_TEST_GOLDENFILES Generate all golden files for test suite
%   This script runs each test category and saves expected outputs
%   to golden files for regression testing.

fprintf('=== Creating Golden Files for Test Suite ===\n\n');

% Add paths
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'src'));
addpath(fullfile(script_dir, '..', 'src', 'autogen'));
addpath(fullfile(script_dir, 'utils'));

% Create goldenfiles directory if it doesn't exist
golden_dir = fullfile(script_dir, 'goldenfiles');
if ~exist(golden_dir, 'dir')
    mkdir(golden_dir);
end

% Metadata for all golden files
metadata = struct();
metadata.version = '1.0';
metadata.date = datestr(now);
metadata.matlab_version = version;

%% 1. Moment Conversions Golden File
fprintf('Creating golden file for moment conversions...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test 1: M2CS4_35 with Gaussian
rho = 1.0; u = 0.1; v = 0.2; w = 0.3; T = 1.5;
M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
[C4, S4] = M2CS4_35(M);
golden_data.tests.gaussian.input = struct('rho', rho, 'u', u, 'v', v, 'w', w, 'T', T);
golden_data.tests.gaussian.output = struct('M', M, 'C4', C4, 'S4', S4);

% Test 2: Correlated case
C200 = 1.2; C020 = 1.5; C002 = 1.8;
C110 = 0.3; C101 = 0.2; C011 = 0.4;
M = InitializeM4_35(rho, u, v, w, C200, C110, C101, C020, C011, C002);
[C4, S4] = M2CS4_35(M);
golden_data.tests.correlated.input = struct('rho', rho, 'u', u, 'v', v, 'w', w, ...
                                             'C200', C200, 'C020', C020, 'C002', C002, ...
                                             'C110', C110, 'C101', C101, 'C011', C011);
golden_data.tests.correlated.output = struct('M', M, 'C4', C4, 'S4', S4);

save(fullfile(golden_dir, 'test_moment_conversions_golden.mat'), 'golden_data');
fprintf('  Saved: test_moment_conversions_golden.mat\n');

%% 2. Initialization Golden File
fprintf('Creating golden file for initialization...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: Various initialization cases
test_cases = {
    struct('rho', 1.0, 'u', 0, 'v', 0, 'w', 0, 'T', 1.0, 'corr', [0,0,0]);
    struct('rho', 1.5, 'u', 0.1, 'v', 0.2, 'w', 0.3, 'T', 1.2, 'corr', [0,0,0]);
    struct('rho', 2.0, 'u', 0, 'v', 0, 'w', 0, 'T', 1.0, 'corr', [0.3,0.2,0.4]);
};

for i = 1:length(test_cases)
    tc = test_cases{i};
    if all(tc.corr == 0)
        M = InitializeM4_35(tc.rho, tc.u, tc.v, tc.w, tc.T, 0, 0, tc.T, 0, tc.T);
    else
        C200 = tc.T; C020 = tc.T; C002 = tc.T;
        C110 = tc.corr(1) * sqrt(C200*C020);
        C101 = tc.corr(2) * sqrt(C200*C002);
        C011 = tc.corr(3) * sqrt(C020*C002);
        M = InitializeM4_35(tc.rho, tc.u, tc.v, tc.w, C200, C110, C101, C020, C011, C002);
    end
    golden_data.tests.(sprintf('case%d', i)).input = tc;
    golden_data.tests.(sprintf('case%d', i)).output = M;
end

save(fullfile(golden_dir, 'test_initialization_golden.mat'), 'golden_data');
fprintf('  Saved: test_initialization_golden.mat\n');

%% 3. Realizability Golden File
fprintf('Creating golden file for realizability...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: S2 realizability
S110_vals = [0.5, 0.99, -0.5];
S101_vals = [0.3, 0.5, 0.4];
S011_vals = [0.4, 0.6, -0.3];

for i = 1:length(S110_vals)
    [S110r, S101r, S011r, S2r] = realizability('S2', S110_vals(i), S101_vals(i), S011_vals(i));
    golden_data.tests.(sprintf('S2_case%d', i)).input = struct('S110', S110_vals(i), ...
                                                                 'S101', S101_vals(i), ...
                                                                 'S011', S011_vals(i));
    golden_data.tests.(sprintf('S2_case%d', i)).output = struct('S110r', S110r, ...
                                                                  'S101r', S101r, ...
                                                                  'S011r', S011r, ...
                                                                  'S2r', S2r);
end

save(fullfile(golden_dir, 'test_realizability_golden.mat'), 'golden_data');
fprintf('  Saved: test_realizability_golden.mat\n');

%% 4. Closures Golden File
fprintf('Creating golden file for closures...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: hyqmom_3D with Gaussian
rho = 1.0; u = 0; v = 0; w = 0; T = 1.0;
M4 = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
[M5, C5, S5] = Moments5_3D(M4);
golden_data.tests.gaussian.input = M4;
golden_data.tests.gaussian.output = struct('M5', M5, 'C5', C5, 'S5', S5);

% Test: closure_and_eigenvalues
mom = [1.0, 0.0, 1.0, 0.0, 3.0];
[Mp, vpmin, vpmax] = closure_and_eigenvalues(mom);
golden_data.tests.closure_1d.input = mom;
golden_data.tests.closure_1d.output = struct('Mp', Mp, 'vpmin', vpmin, 'vpmax', vpmax);

save(fullfile(golden_dir, 'test_closures_golden.mat'), 'golden_data');
fprintf('  Saved: test_closures_golden.mat\n');

%% 5. Flux and Eigenvalues Golden File
fprintf('Creating golden file for flux and eigenvalues...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: Flux computation
rho = 1.0; u = 0.1; v = 0.2; w = 0.3; T = 1.0;
M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
[Fx, Fy, Fz, M_real] = Flux_closure35_and_realizable_3D(M, 1, 0.5);
golden_data.tests.flux_basic.input = M;
golden_data.tests.flux_basic.output = struct('Fx', Fx, 'Fy', Fy, 'Fz', Fz, 'M_real', M_real);

% Test: Eigenvalues
[v6min_x, v6max_x, ~] = eigenvalues6_hyperbolic_3D(M, 1, 1, 0.5);
[v6min_y, v6max_y, ~] = eigenvalues6_hyperbolic_3D(M, 2, 1, 0.5);
golden_data.tests.eigenvalues.input = M;
golden_data.tests.eigenvalues.output = struct('v6min_x', v6min_x, 'v6max_x', v6max_x, ...
                                                'v6min_y', v6min_y, 'v6max_y', v6max_y);

save(fullfile(golden_dir, 'test_flux_eigenvalues_golden.mat'), 'golden_data');
fprintf('  Saved: test_flux_eigenvalues_golden.mat\n');

%% 6. Numerical Schemes Golden File
fprintf('Creating golden file for numerical schemes...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: pas_HLL
Np = 20; Nmom = 35;
M = ones(Np, Nmom);
F = ones(Np, Nmom);
dt = 0.01; dx = 0.1;
vpmin = -ones(Np, 1); vpmax = ones(Np, 1);
Mp = pas_HLL(M, F, dt, dx, vpmin, vpmax, true, true);
golden_data.tests.pas_hll.input = struct('M', M, 'F', F, 'dt', dt, 'dx', dx, ...
                                          'vpmin', vpmin, 'vpmax', vpmax);
golden_data.tests.pas_hll.output = Mp;

save(fullfile(golden_dir, 'test_numerical_schemes_golden.mat'), 'golden_data');
fprintf('  Saved: test_numerical_schemes_golden.mat\n');

%% 7. Collision Golden File
fprintf('Creating golden file for collision...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: collision35
rho = 1.5; u = 0.1; v = 0.2; w = 0.3; T = 1.0;
M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
dt = 0.01; Kn = 1.0;
M_out = collision35(M, dt, Kn);
golden_data.tests.basic.input = struct('M', M, 'dt', dt, 'Kn', Kn);
golden_data.tests.basic.output = M_out;

save(fullfile(golden_dir, 'test_collision_golden.mat'), 'golden_data');
fprintf('  Saved: test_collision_golden.mat\n');

%% 8. Diagnostics Golden File
fprintf('Creating golden file for diagnostics...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: check2D
S30 = 0.5; S40 = 3.5; S11 = 0.3;
S21 = 0.2; S31 = 0.1; S12 = 0.15;
S22 = 1.2; S03 = 0.4; S13 = 0.1; S04 = 3.2;
[S30_out, S40_out, S11_out, S21_out, S31_out, S12_out, S22_out, S03_out, S13_out, S04_out] = ...
    diagnostics('check2D', S30, S40, S11, S21, S31, S12, S22, S03, S13, S04);
golden_data.tests.check2d.input = struct('S30', S30, 'S40', S40, 'S11', S11, ...
                                          'S21', S21, 'S31', S31, 'S12', S12, ...
                                          'S22', S22, 'S03', S03, 'S13', S13, 'S04', S04);
golden_data.tests.check2d.output = struct('S30', S30_out, 'S40', S40_out, 'S11', S11_out, ...
                                           'S21', S21_out, 'S31', S31_out, 'S12', S12_out, ...
                                           'S22', S22_out, 'S03', S03_out, 'S13', S13_out, 'S04', S04_out);

save(fullfile(golden_dir, 'test_diagnostics_golden.mat'), 'golden_data');
fprintf('  Saved: test_diagnostics_golden.mat\n');

%% 9. Grid Utils Golden File
fprintf('Creating golden file for grid utils...\n');
golden_data = struct();
golden_data.metadata = metadata;
golden_data.tests = struct();

% Test: grid setup
Np = 20;
grid = grid_utils('setup', Np, -0.5, 0.5, -0.5, 0.5);
golden_data.tests.grid_setup.input = struct('Np', Np);
golden_data.tests.grid_setup.output = grid;

save(fullfile(golden_dir, 'test_grid_utils_golden.mat'), 'golden_data');
fprintf('  Saved: test_grid_utils_golden.mat\n');

