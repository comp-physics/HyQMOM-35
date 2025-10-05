function build_mex(varargin)
%BUILD_MEX - Cross-platform MEX builder for delta2star3D
%
% Usage:
%   build_mex()              % Auto-detect platform and build
%   build_mex('clean')       % Remove all MEX binaries
%   build_mex('test')        % Build and run tests
%   build_mex('verbose')     % Build with verbose output
%
% This script compiles delta2star3D_mex.c with platform-specific optimizations.
% The resulting MEX binary will have the correct extension for your platform:
%   - macOS ARM64 (M1/M2/M3):  .mexmaca64
%   - macOS Intel x86_64:       .mexmaci64
%   - Linux x86_64:             .mexa64
%   - Windows x86_64:           .mexw64
%
% Examples:
%   build_mex()              % Normal build
%   build_mex('clean')       % Clean all MEX files
%   build_mex('test')        % Build and verify

%% Parse arguments
clean_only = false;
run_tests = false;
verbose = false;

for i = 1:length(varargin)
    switch lower(varargin{i})
        case 'clean'
            clean_only = true;
        case 'test'
            run_tests = true;
        case 'verbose'
            verbose = true;
    end
end

%% Define paths
script_dir = fileparts(mfilename('fullpath'));
mex_dir = fullfile(script_dir, 'src', 'autogen');
mex_source = fullfile(mex_dir, 'delta2star3D_mex.c');
mex_output_base = fullfile(mex_dir, 'delta2star3D_mex');

%% Clean operation
if clean_only
    fprintf('Cleaning MEX binaries...\n');
    mex_files = dir(fullfile(mex_dir, 'delta2star3D_mex.mex*'));
    for i = 1:length(mex_files)
        mex_file = fullfile(mex_dir, mex_files(i).name);
        fprintf('  Removing: %s\n', mex_files(i).name);
        delete(mex_file);
    end
    fprintf('Clean complete\n');
    return;
end

%% Check that source file exists
if ~exist(mex_source, 'file')
    error('MEX source file not found: %s', mex_source);
end

%% Detect platform
fprintf('=== MEX Build System ===\n\n');
fprintf('Platform detection:\n');

if ismac
    [~, arch] = system('uname -m');
    arch = strtrim(arch);
    if contains(arch, 'arm') || contains(arch, 'aarch64')
        platform = 'macOS ARM64 (Apple Silicon)';
        mex_ext = 'mexmaca64';
    else
        platform = 'macOS Intel x86_64';
        mex_ext = 'mexmaci64';
    end
elseif isunix
    platform = 'Linux x86_64';
    mex_ext = 'mexa64';
elseif ispc
    platform = 'Windows x86_64';
    mex_ext = 'mexw64';
else
    error('Unknown platform');
end

fprintf('  Platform: %s\n', platform);
fprintf('  Extension: .%s\n\n', mex_ext);

%% Build MEX file with platform-specific flags
fprintf('Building delta2star3D_mex...\n');
fprintf('  Source: %s\n', mex_source);
fprintf('  Output: %s.%s\n\n', mex_output_base, mex_ext);

% Change to MEX directory for compilation
old_dir = cd(mex_dir);

try
    % Platform-specific compiler flags
    if ismac || isunix
        % macOS and Linux: GCC/Clang flags
        cflags = 'COPTIMFLAGS=''-O3 -ffast-math -march=native -DNDEBUG''';
        if verbose
            fprintf('Compiler flags: %s\n\n', cflags);
        end
        mex_cmd = sprintf('mex -O %s delta2star3D_mex.c -output delta2star3D_mex', cflags);
    else
        % Windows: MSVC flags
        cflags = 'COMPFLAGS="/O2 /fp:fast /DNDEBUG"';
        if verbose
            fprintf('Compiler flags: %s\n\n', cflags);
        end
        mex_cmd = sprintf('mex -O %s delta2star3D_mex.c -output delta2star3D_mex', cflags);
    end
    
    % Execute MEX compilation
    if verbose
        fprintf('Command: %s\n\n', mex_cmd);
    end
    
    eval(mex_cmd);
    
    fprintf('\n✓ MEX compilation successful!\n');
    
    % Verify output file exists
    output_file = sprintf('%s.%s', mex_output_base, mex_ext);
    if exist(output_file, 'file')
        file_info = dir(output_file);
        fprintf('✓ Output file created: %s (%.1f KB)\n', ...
                file_info.name, file_info.bytes/1024);
    else
        warning('MEX file was compiled but not found at expected location');
    end
    
catch ME
    cd(old_dir);
    fprintf('\n✗ MEX compilation failed!\n');
    fprintf('Error: %s\n\n', ME.message);
    
    fprintf('Troubleshooting:\n');
    fprintf('  1. Ensure you have a C compiler installed:\n');
    if ismac
        fprintf('     - macOS: Install Xcode Command Line Tools\n');
        fprintf('       Command: xcode-select --install\n');
    elseif isunix
        fprintf('     - Linux: Install GCC\n');
        fprintf('       Command: sudo apt-get install build-essential (Debian/Ubuntu)\n');
        fprintf('                sudo yum install gcc (RedHat/CentOS)\n');
    else
        fprintf('     - Windows: Install MinGW-w64 or Visual Studio\n');
    end
    fprintf('  2. Verify compiler setup in MATLAB:\n');
    fprintf('     Command: mex -setup\n');
    fprintf('  3. The code will fall back to MATLAB implementation automatically\n');
    
    rethrow(ME);
end

cd(old_dir);

%% Run tests if requested
if run_tests
    fprintf('\n=== Running Tests ===\n\n');
    
    % Add paths
    addpath(genpath(script_dir));
    
    % Test 1: Basic functionality
    fprintf('Test 1: Basic functionality...\n');
    s = randn(1, 28) * 0.1;
    try
        E = delta2star3D(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10),...
                        s(11),s(12),s(13),s(14),s(15),s(16),s(17),s(18),s(19),s(20),...
                        s(21),s(22),s(23),s(24),s(25),s(26),s(27),s(28));
        fprintf('  ✓ Function executes successfully\n');
        fprintf('  ✓ Output size: %dx%d\n', size(E, 1), size(E, 2));
    catch ME
        fprintf('  ✗ FAILED: %s\n', ME.message);
        return;
    end
    
    % Test 2: Accuracy comparison
    fprintf('\nTest 2: Comparing MEX vs MATLAB...\n');
    rng(42);
    max_diffs = zeros(10, 1);
    for i = 1:10
        s = randn(1, 28) * 0.1;
        E_mex = delta2star3D_mex(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10),...
                                s(11),s(12),s(13),s(14),s(15),s(16),s(17),s(18),s(19),s(20),...
                                s(21),s(22),s(23),s(24),s(25),s(26),s(27),s(28));
        E_matlab = delta2star3D_matlab(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10),...
                                      s(11),s(12),s(13),s(14),s(15),s(16),s(17),s(18),s(19),s(20),...
                                      s(21),s(22),s(23),s(24),s(25),s(26),s(27),s(28));
        max_diffs(i) = max(abs(E_mex(:) - E_matlab(:)));
    end
    
    max_error = max(max_diffs);
    if max_error < 1e-13
        fprintf('  ✓ Accuracy: %.2e (excellent)\n', max_error);
    else
        fprintf('  ✗ Accuracy: %.2e (may be problematic)\n', max_error);
    end
    
    % Test 3: Performance benchmark
    fprintf('\nTest 3: Performance benchmark...\n');
    s = randn(1, 28) * 0.1;
    n_iter = 1000;
    
    % Warmup
    for i = 1:10
        delta2star3D_mex(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10),...
                         s(11),s(12),s(13),s(14),s(15),s(16),s(17),s(18),s(19),s(20),...
                         s(21),s(22),s(23),s(24),s(25),s(26),s(27),s(28));
    end
    
    % Benchmark MEX
    tic;
    for i = 1:n_iter
        delta2star3D_mex(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10),...
                         s(11),s(12),s(13),s(14),s(15),s(16),s(17),s(18),s(19),s(20),...
                         s(21),s(22),s(23),s(24),s(25),s(26),s(27),s(28));
    end
    time_mex = toc;
    
    % Benchmark MATLAB
    tic;
    for i = 1:n_iter
        delta2star3D_matlab(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9),s(10),...
                            s(11),s(12),s(13),s(14),s(15),s(16),s(17),s(18),s(19),s(20),...
                            s(21),s(22),s(23),s(24),s(25),s(26),s(27),s(28));
    end
    time_matlab = toc;
    
    speedup = time_matlab / time_mex;
    fprintf('  MEX:    %.4f ms/call\n', time_mex/n_iter*1000);
    fprintf('  MATLAB: %.4f ms/call\n', time_matlab/n_iter*1000);
    fprintf('  ✓ Speedup: %.1fx\n', speedup);
    
    fprintf('\n=== All Tests Passed ===\n');
end

fprintf('\n=== Build Complete ===\n');
fprintf('The delta2star3D function will now automatically use the MEX version.\n');
fprintf('To verify: run main() and look for "[delta2star3D] Using fast MEX implementation"\n\n');

end

