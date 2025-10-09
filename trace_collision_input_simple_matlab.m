% Ultra-simple trace: Just run simulation and save moments before/after collision

clear all;
close all;

setup_paths;

fprintf('===== SIMPLE TRACE - MATLAB =====\n\n');

% Instrument collision35 to save its input
global COLLISION_INPUT;
COLLISION_INPUT = [];

% Save original
orig_collision = which('collision35');

% Create temporary instrumented version
temp_dir = tempname;
mkdir(temp_dir);
addpath(temp_dir, '-begin');

% Write instrumented version
fid = fopen(fullfile(temp_dir, 'collision35.m'), 'w');
fprintf(fid, 'function M_out = collision35(M_in, dt, Kn)\n');
fprintf(fid, '  global COLLISION_INPUT;\n');
fprintf(fid, '  if isempty(COLLISION_INPUT)\n');
fprintf(fid, '    COLLISION_INPUT = M_in;\n');
fprintf(fid, '    fprintf(''Captured collision input: M210=%%e, M130=%%e\\n'', M_in(8), M_in(14));\n');
fprintf(fid, '  end\n');
fprintf(fid, '  % Call original\n');
fprintf(fid, '  load(''%s'');\n', orig_collision);
fprintf(fid, '  M_out = collision35_original(M_in, dt, Kn);\n');
fprintf(fid, 'end\n');
fclose(fid);

% Run simulation for one step
params = struct('Np', 20, 'Kn', 1.0, 'Ma', 0.0, 'tmax', 0.05, 'flag2D', 0, 'CFL', 0.5);

fprintf('Running simulation...\n');
result = simulation_runner(params);

% Check what collision saw
if ~isempty(COLLISION_INPUT)
    fprintf('\n=== COLLISION INPUT (first call) ===\n');
    fprintf('M210 = %.15e\n', COLLISION_INPUT(8));
    fprintf('M130 = %.15e\n', COLLISION_INPUT(14));
else
    fprintf('\nNo collision input captured\n');
end

% Cleanup
rmpath(temp_dir);
rmdir(temp_dir, 's');

fprintf('\n===== END MATLAB TRACE =====\n');
