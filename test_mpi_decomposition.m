function test_mpi_decomposition()
%TEST_MPI_DECOMPOSITION Comprehensive test suite for MPI domain decomposition
%   This script tests the MPI implementation with various configurations
%   and verifies correctness against serial reference solutions.

    % Add src directory to path
    script_dir = fileparts(mfilename('fullpath'));
    src_dir = fullfile(script_dir, 'src');
    if exist(src_dir, 'dir')
        addpath(src_dir);
    end

    fprintf('\n╔════════════════════════════════════════════════════╗\n');
    fprintf('║  MPI Domain Decomposition Test Suite              ║\n');
    fprintf('╚════════════════════════════════════════════════════╝\n\n');

    % Test configurations
    configs = {
        struct('name', 'Small grid, 4 workers',  'np', 16, 'workers', 4,  'steps', 3, 'nv', 2)
        struct('name', 'Medium grid, 4 workers', 'np', 32, 'workers', 4,  'steps', 5, 'nv', 2)
        struct('name', 'Large grid, 4 workers',  'np', 64, 'workers', 4,  'steps', 5, 'nv', 3)
        struct('name', 'Medium grid, 9 workers', 'np', 36, 'workers', 9,  'steps', 5, 'nv', 2)
        struct('name', 'Small grid, 1 worker',   'np', 16, 'workers', 1,  'steps', 3, 'nv', 2)
    };

    results = cell(length(configs), 1);
    
    for i = 1:length(configs)
        cfg = configs{i};
        fprintf('Test %d/%d: %s\n', i, length(configs), cfg.name);
        fprintf('  Grid: %d×%d, Workers: %d, Steps: %d, Variables: %d\n', ...
                cfg.np, cfg.np, cfg.workers, cfg.steps, cfg.nv);
        
        try
            err = run_single_test(cfg.np, cfg.steps, cfg.nv, cfg.workers);
            results{i} = struct('config', cfg, 'error', err, 'passed', err < 1e-12);
            
            if err < 1e-12
                fprintf('  ✓ PASS - Error: %.3e\n\n', err);
            else
                fprintf('  ✗ FAIL - Error: %.3e (exceeds tolerance)\n\n', err);
            end
        catch ME
            fprintf('  ✗ ERROR: %s\n\n', ME.message);
            results{i} = struct('config', cfg, 'error', inf, 'passed', false, 'exception', ME);
        end
    end
    
    % Summary
    fprintf('╔════════════════════════════════════════════════════╗\n');
    fprintf('║  Test Summary                                      ║\n');
    fprintf('╚════════════════════════════════════════════════════╝\n\n');
    
    n_passed = sum(cellfun(@(r) r.passed, results));
    n_total = length(results);
    
    fprintf('Passed: %d / %d tests\n', n_passed, n_total);
    
    if n_passed == n_total
        fprintf('\n🎉 All tests PASSED!\n\n');
    else
        fprintf('\n⚠️  Some tests FAILED. Check output above for details.\n\n');
        
        % Show failed tests
        fprintf('Failed tests:\n');
        for i = 1:length(results)
            if ~results{i}.passed
                fprintf('  - %s (error: %.3e)\n', results{i}.config.name, results{i}.error);
            end
        end
        fprintf('\n');
    end
end

function err = run_single_test(np, steps, nv, num_workers)
    % Start/restart parallel pool
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', num_workers);
    elseif pool.NumWorkers ~= num_workers
        delete(pool);
        parpool('local', num_workers);
    end

    halo = 1;
    bc = struct('type', 'copy');

    globalOut = [];
    baseField = [];
    
    spmd
        decomp = setup_mpi_cartesian_2d(np, halo);
        
        if labindex == 1
            baseField = zeros(np, np, nv);
            for v = 1:nv
                [X, Y] = ndgrid(1:np, 1:np);
                % Use different pattern for variety
                baseField(:,:,v) = sin(2*pi*X/np).*cos(2*pi*Y/np) + ...
                                   0.5*cos(3*pi*X/np).*sin(3*pi*Y/np) + v;
            end
        end
        baseField = labBroadcast(1, baseField);

        nx = decomp.local_size(1);
        ny = decomp.local_size(2);
        A  = zeros(nx+2*halo, ny+2*halo, nv);
        Ap = A;

        i0i1 = decomp.istart_iend;
        j0j1 = decomp.jstart_jend;
        A(halo+1:halo+nx, halo+1:halo+ny, :) = baseField(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :);

        for s = 1:steps
            A = halo_exchange_2d(A, decomp, bc);

            for v = 1:nv
                C = A(:,:,v);
                Ap(:,:,v) = C;
                Ap(halo+1:halo+nx, halo+1:halo+ny, v) = 0.25 * ( ...
                    C(halo+2:halo+nx+1, halo+1:halo+ny) + ...
                    C(halo:halo+nx-1,   halo+1:halo+ny) + ...
                    C(halo+1:halo+nx,   halo+2:halo+ny+1) + ...
                    C(halo+1:halo+nx,   halo:halo+ny-1) ...
                );
            end
            A = Ap;
        end

        gathered = labSendReceive(1, 1, A(halo+1:halo+nx, halo+1:halo+ny, :));
        if labindex == 1
            globalOut = zeros(np, np, nv);
            globalOut(decomp.istart_iend(1):decomp.istart_iend(2), ...
                      decomp.jstart_jend(1):decomp.jstart_jend(2), :) = gathered;
            for src = 2:numlabs
                blk = labReceive(src);
                Px = decomp.dims(1);
                r  = src - 1;
                rx = mod(r, Px);
                ry = floor(r / Px);
                [~, i0_s, i1_s] = local_block(np, decomp.dims(1), rx);
                [~, j0_s, j1_s] = local_block(np, decomp.dims(2), ry);
                globalOut(i0_s:i1_s, j0_s:j1_s, :) = blk;
            end
        else
            labSend(A(halo+1:halo+nx, halo+1:halo+ny, :), 1);
        end
    end

    if isa(globalOut, 'Composite')
        G = globalOut{1};
    else
        G = globalOut;
    end

    % Serial reference
    Ref = baseField;
    for s = 1:steps
        R = Ref;
        for v = 1:nv
            C = Ref(:,:,v);
            C([1,np],:) = C([2,np-1],:);
            C(:,[1,np]) = C(:,[2,np-1]);
            R(:,:,v) = 0.25 * ( ...
                C([2:np, np], :) + C([1, 1:np-1], :) + ...
                C(:, [2:np, np]) + C(:, [1, 1:np-1]) ...
            );
        end
        Ref = R;
    end

    err = max(abs(G(:) - Ref(:)));
end

function [n_local, i0, i1] = local_block(n, P, r)
    base = floor(n/P);
    remn = mod(n, P);
    if r < remn
        n_local = base + 1;
        i0 = r*(base+1) + 1;
    else
        n_local = base;
        i0 = remn*(base+1) + (r-remn)*base + 1;
    end
    i1 = i0 + n_local - 1;
end

