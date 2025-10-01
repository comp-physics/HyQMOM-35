function example_mpi_solver_loop(np, steps, nv, num_workers)
%EXAMPLE_MPI_SOLVER_LOOP Minimal example using 2D domain decomposition + halo exchanges.
%   example_mpi_solver_loop(np, steps, nv, num_workers)
%
% This uses a toy 5-point stencil to demonstrate correctness of MPI exchanges.
%
% Inputs:
%   np          - global grid size (default: 32)
%   steps       - number of stencil iterations (default: 5)
%   nv          - number of variables (default: 2)
%   num_workers - number of parallel workers (default: 4)
%
% Example:
%   example_mpi_solver_loop()           % Run with defaults
%   example_mpi_solver_loop(64, 10, 3, 4)  % Custom parameters

    if nargin < 1, np = 32; end
    if nargin < 2, steps = 5; end
    if nargin < 3, nv = 2; end
    if nargin < 4, num_workers = 4; end

    % Add src directory to path
    script_dir = fileparts(mfilename('fullpath'));
    src_dir = fullfile(script_dir, 'src');
    if exist(src_dir, 'dir')
        addpath(src_dir);
    end

    % Start parallel pool if needed
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('Starting parallel pool with %d workers...\n', num_workers);
        parpool('local', num_workers);
    elseif pool.NumWorkers ~= num_workers
        delete(pool);
        fprintf('Restarting parallel pool with %d workers...\n', num_workers);
        parpool('local', num_workers);
    end

    halo = 1;
    bc = struct('type', 'copy');

    fprintf('\n=== MPI Domain Decomposition Example ===\n');
    fprintf('Grid size: %d x %d\n', np, np);
    fprintf('Workers: %d\n', num_workers);
    fprintf('Stencil iterations: %d\n', steps);
    fprintf('Variables: %d\n\n', nv);

    % Initialize a known global field on lab 1, broadcast to all labs
    globalOut = [];
    baseField = [];
    
    spmd
        decomp = setup_mpi_cartesian_2d(np, halo);
        
        if labindex == 1
            fprintf('Process grid: %d x %d\n', decomp.dims(1), decomp.dims(2));
            fprintf('Local subdomain size: %d x %d (per worker)\n\n', ...
                    decomp.local_size(1), decomp.local_size(2));
            baseField = zeros(np, np, nv);
            for v = 1:nv
                [X, Y] = ndgrid(1:np, 1:np);
                baseField(:,:,v) = sin(2*pi*X/np).*cos(2*pi*Y/np) + v;
            end
        end
        baseField = labBroadcast(1, baseField);

        % Allocate local with halos
        nx = decomp.local_size(1);
        ny = decomp.local_size(2);
        A  = zeros(nx+2*halo, ny+2*halo, nv);
        Ap = A;

        % Scatter interior slice from global field
        i0i1 = decomp.istart_iend;
        j0j1 = decomp.jstart_jend;
        A(halo+1:halo+nx, halo+1:halo+ny, :) = baseField(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :);

        % Time stepping with 5-point stencil
        for s = 1:steps
            % Exchange halos
            A = halo_exchange_2d(A, decomp, bc);

            % Toy 5-point stencil update on interior only
            for v = 1:nv
                C = A(:,:,v);
                Ap(:,:,v) = C; % copy
                Ap(halo+1:halo+nx, halo+1:halo+ny, v) = 0.25 * ( ...
                    C(halo+2:halo+nx+1, halo+1:halo+ny) + ...   % i+1
                    C(halo:halo+nx-1,   halo+1:halo+ny) + ...   % i-1
                    C(halo+1:halo+nx,   halo+2:halo+ny+1) + ... % j+1
                    C(halo+1:halo+nx,   halo:halo+ny-1) ...     % j-1
                );
            end
            A = Ap;
        end

        % Gather local interiors back to lab 1 for verification
        gathered = labSendReceive(1, 1, A(halo+1:halo+nx, halo+1:halo+ny, :));
        if labindex == 1
            globalOut = zeros(np, np, nv);
            % Place own block first
            globalOut(decomp.istart_iend(1):decomp.istart_iend(2), ...
                      decomp.jstart_jend(1):decomp.jstart_jend(2), :) = gathered;
            % Receive from others
            for src = 2:numlabs
                blk = labReceive(src);
                % Compute src decomp (deterministic mapping)
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

    % Reference serial result to verify (same toy stencil)
    fprintf('Computing serial reference solution...\n');
    Ref = baseField;
    for s = 1:steps
        R = Ref;
        for v = 1:nv
            C = Ref(:,:,v);
            % Apply 'copy' BC
            C([1,np],:) = C([2,np-1],:);
            C(:,[1,np]) = C(:,[2,np-1]);
            % 5-point stencil
            R(:,:,v) = 0.25 * ( ...
                C([2:np, np], :) + C([1, 1:np-1], :) + ...
                C(:, [2:np, np]) + C(:, [1, 1:np-1]) ...
            );
        end
        Ref = R;
    end

    err = max(abs(G(:) - Ref(:)));
    fprintf('\n=== VERIFICATION RESULTS ===\n');
    fprintf('Max absolute error: %.3e\n', err);
    if err < 1e-12
        fprintf('✓ PASS: MPI domain decomposition working correctly!\n\n');
    else
        fprintf('✗ FAIL: Error exceeds tolerance\n\n');
    end
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

