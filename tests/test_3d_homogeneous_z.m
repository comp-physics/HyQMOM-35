function tests = test_3d_homogeneous_z
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(script_dir);  % Add parent directory for main.m
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-10;  % Tight tolerance for identical physics
    testCase.TestData.has_pct = license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
end

function test_3d_basic_run(testCase)
    % Basic test that 3D simulation runs without errors
    % Uses inhomogeneous-in-z IC (jets in lower half of z domain)
    
    if ~testCase.TestData.has_pct
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required');
        return;
    end
    
    fprintf('\n=== Test: 3D Basic Run ===\n');
    
    % Small grid for fast testing
    Np = 20;  % Needs to be >= 10*num_workers for MPI decomposition
    Nz = 4;
    tmax = 0.001;  % Very short time
    num_workers = 1;  % Single worker for faster testing
    
    fprintf('Running 3D simulation (Np=%d, Nz=%d)...\n', Np, Nz);
    results_3d = main(Np, tmax, false, num_workers, false, 1, false, Nz);
    
    % Basic sanity checks
    verifyTrue(testCase, results_3d.parameters.time_steps >= 1, 'Should take at least 1 time step');
    verifyEqual(testCase, size(results_3d.moments.M, 1), Np, 'X grid size should match');
    verifyEqual(testCase, size(results_3d.moments.M, 2), Np, 'Y grid size should match');
    verifyEqual(testCase, size(results_3d.moments.M, 3), Nz, 'Z grid size should match');
    
    % Check that density remains positive
    M_3d = results_3d.moments.M;
    verifyTrue(testCase, all(M_3d(:,:,:,1) > 0, 'all'), 'Density should remain positive');
    
    % Check that z-variation exists (jets only in lower z)
    M_lower_z = M_3d(:,:,1,1);  % Density at z-slice 1 (lower, should have jets)
    M_upper_z = M_3d(:,:,end,1);  % Density at z-slice end (upper, should be background)
    mean_lower = mean(M_lower_z(:));
    mean_upper = mean(M_upper_z(:));
    verifyGreaterThan(testCase, mean_lower, mean_upper, ...
                      'Lower z should have higher mean density (jets present)');
    
    fprintf('PASS: 3D simulation runs correctly with z-variation\n');
end

function test_z_slice_in_jet_region_uniform(testCase)
    % Test that z-slices within the jet region (z < 0) are similar to each other
    % This validates that the physics is working correctly in z
    
    if ~testCase.TestData.has_pct
        fprintf('SKIPPED: Parallel Computing Toolbox not available\n');
        assumeFail(testCase, 'Parallel Computing Toolbox required');
        return;
    end
    
    fprintf('\n=== Test: Z-slices in jet region should be similar ===\n');
    
    Np = 20;
    Nz = 8;  % More z-slices to see gradient
    tmax = 0.001;
    num_workers = 1;  % Single worker for faster testing
    
    fprintf('Running 3D simulation with Nz=%d...\n', Nz);
    results_3d = main(Np, tmax, false, num_workers, false, 1, false, Nz);
    
    M_3d = results_3d.moments.M;  % (Np, Np, Nz, Nmom)
    
    % Find z-coordinates
    zm = results_3d.grid.zm;
    
    % Find z-slices in jet region (z < 0)
    jet_slices = find(zm < 0);
    background_slices = find(zm >= 0);
    
    fprintf('Jet region: z-slices %s\n', mat2str(jet_slices));
    fprintf('Background region: z-slices %s\n', mat2str(background_slices));
    
    if length(jet_slices) >= 2
        % Compare jet slices - they should be nearly identical (homogeneous in z within jet region)
        k1 = jet_slices(1);
        k2 = jet_slices(end);
        diff_jet = M_3d(:,:,k1,:) - M_3d(:,:,k2,:);
        max_diff_jet = max(abs(diff_jet(:)));
        rel_diff_jet = max_diff_jet / (max(abs(M_3d(:,:,k1,1)), [], 'all') + eps);
        fprintf('Jet region z-slices %d vs %d: max diff = %.3e, relative = %.3e\n', ...
                k1, k2, max_diff_jet, rel_diff_jet);
        
        % Should be reasonably similar (allow for some z-transport effects)
        % Note: Even though there's no MPI decomposition in z, physics can cause z-variation
        verifyLessThan(testCase, rel_diff_jet, 0.5, ...
                       'Jet region z-slices should be reasonably similar (< 50% relative difference)');
    end
    
    if length(background_slices) >= 2
        % Compare background slices
        k1 = background_slices(1);
        k2 = background_slices(end);
        diff_bg = M_3d(:,:,k1,:) - M_3d(:,:,k2,:);
        max_diff_bg = max(abs(diff_bg(:)));
        rel_diff_bg = max_diff_bg / (max(abs(M_3d(:,:,k1,1)), [], 'all') + eps);
        fprintf('Background region z-slices %d vs %d: max diff = %.3e, relative = %.3e\n', ...
                k1, k2, max_diff_bg, rel_diff_bg);
        
        % Should be reasonably similar (background can evolve due to z-transport)
        verifyLessThan(testCase, rel_diff_bg, 2.0, ...
                       'Background region z-slices should be reasonably similar (< 200% relative difference)');
    end
    
    fprintf('PASS: Z-slices within each region are similar\n');
end


