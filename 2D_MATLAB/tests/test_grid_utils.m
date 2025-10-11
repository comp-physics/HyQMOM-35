function tests = test_grid_utils
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-12;
end

function test_grid_utils_setup(testCase)
    Np = 20;
    xmin = -0.5; xmax = 0.5;
    ymin = -0.5; ymax = 0.5;
    
    grid = grid_utils('setup', Np, xmin, xmax, ymin, ymax);
    
    verifyEqual(testCase, length(grid.x), Np+1, 'x should have Np+1 points');
    verifyEqual(testCase, length(grid.y), Np+1, 'y should have Np+1 points');
    verifyEqual(testCase, length(grid.xm), Np, 'xm should have Np points');
    verifyEqual(testCase, length(grid.ym), Np, 'ym should have Np points');
    
    verifyEqual(testCase, grid.dx, (xmax-xmin)/Np, 'AbsTol', testCase.TestData.tol, 'dx correct');
    verifyEqual(testCase, grid.dy, (ymax-ymin)/Np, 'AbsTol', testCase.TestData.tol, 'dy correct');
    
    fprintf('PASS: grid_utils setup test\n');
end

function test_grid_moment_processor(testCase)
    Np = 10;
    Nmom = 35;
    
    M = zeros(Np, Np, Nmom);
    for i = 1:Np
        for j = 1:Np
            rho = 1.0 + 0.1*sin(i) + 0.1*cos(j);
            M(i,j,:) = InitializeM4_35(rho, 0, 0, 0, 1, 0, 0, 1, 0, 1);
        end
    end
    
    [C, S] = grid_moment_processor(M, @M2CS4_35);
    
    verifyEqual(testCase, size(C), size(M), 'C should have same size as M');
    verifyEqual(testCase, size(S), size(M), 'S should have same size as M');
    verifyTrue(testCase, all(isfinite(C(:))), 'All C values should be finite');
    verifyTrue(testCase, all(isfinite(S(:))), 'All S values should be finite');
    
    fprintf('PASS: grid_moment_processor test\n');
end

function test_moment_idx(testCase)
    idx = moment_idx('M000');
    verifyEqual(testCase, idx, 1, 'M000 should be at index 1');
    
    idx = moment_idx('M100');
    verifyEqual(testCase, idx, 2, 'M100 should be at index 2');
    
    idx = moment_idx('M110');
    verifyEqual(testCase, idx, 7, 'M110 should be at index 7');
    
    fprintf('PASS: moment_idx test\n');
end
