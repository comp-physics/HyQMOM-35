function [M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir)
% Run parallel SPMD simulation loop
% This is the core time-stepping loop extracted from main.m
%
% Inputs:
%   params     - Struct with all simulation parameters
%   script_dir - Path to project root for worker path setup
%
% Outputs:
%   M_final    - Final moment array (global, gathered to rank 1)
%   final_time - Final simulation time
%   time_steps - Number of time steps taken
%   grid_out   - Grid structure

% Setup grid and initial conditions (outside spmd)
grid = grid_utils('setup', params.Np, -0.5, 0.5, -0.5, 0.5);
M_global = grid_utils('crossing_jets_ic', params.Np, params.Nmom, ...
                      params.rhol, params.rhor, params.Ma, params.T, ...
                      params.r110, params.r101, params.r011);

% MPI parallel execution
spmd
    % Workers need src/ directory in their path
    setup_paths(script_dir);
    
    % Define all constants directly inside spmd block
    halo = 2;  % Required for pas_HLL stencil
    bc = struct('type', 'copy');
    
    % Unpack parameters for workers (avoid Composite overhead in loops)
    Np = params.Np;
    tmax = params.tmax;
    Kn_worker = params.Kn;
    Ma_worker = params.Ma;
    flag2D_worker = params.flag2D;
    CFL_worker = params.CFL;
    dx_worker = params.dx;
    dy_worker = params.dy;
    Nmom_worker = params.Nmom;
    Nmom5_worker = params.Nmom5;
    nnmax_worker = params.nnmax;
    dtmax_worker = params.dtmax;
    
    % Setup domain decomposition
    decomp = setup_mpi_cartesian_2d(Np, halo);
    nx = decomp.local_size(1);
    ny = decomp.local_size(2);
    
    % Allocate local arrays with halos
    M = zeros(nx+2*halo, ny+2*halo, Nmom_worker);
    Mnp = M;
    Fx = zeros(nx+2*halo, ny+2*halo, Nmom_worker);
    Fy = zeros(nx+2*halo, ny+2*halo, Nmom_worker);
    
    % Wave speed arrays (interior only, no halos needed)
    vpxmin = zeros(nx, ny);
    vpxmax = zeros(nx, ny);
    vpymin = zeros(nx, ny);
    vpymax = zeros(nx, ny);
    v5xmin = zeros(nx, ny);
    v5xmax = zeros(nx, ny);
    v5ymin = zeros(nx, ny);
    v5ymax = zeros(nx, ny);
    v6xmin = zeros(nx, ny);
    v6xmax = zeros(nx, ny);
    v6ymin = zeros(nx, ny);
    v6ymax = zeros(nx, ny);
    
    % Scatter initial conditions to local subdomains
    i0i1 = decomp.istart_iend;
    j0j1 = decomp.jstart_jend;
    M(halo+1:halo+nx, halo+1:halo+ny, :) = M_global(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :);
    
    % Initial halo exchange
    M = halo_exchange_2d(M, decomp, bc);
    
    % Time evolution
    t = 0.0;
    nn = 0;
    
    while t < tmax && nn < nnmax_worker
        nn = nn + 1;
        step_start_time = tic;  % Start timing this timestep
        
        % Compute fluxes and wave speeds for interior cells
        for i = 1:nx
            for j = 1:ny
                % Access with halo offset
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D_worker, Ma_worker);
                [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D_worker, Ma_worker);
                [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D_worker, Ma_worker);
                
                Fx(ih, jh, :) = Mx;
                Fy(ih, jh, :) = My;
                Mnp(ih, jh, :) = Mr;
                
                [~, v5xmin(i,j), v5xmax(i,j)] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
                [~, v5ymin(i,j), v5ymax(i,j)] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
                
                vpxmin(i,j) = min(v5xmin(i,j), v6xmin(i,j));
                vpxmax(i,j) = max(v5xmax(i,j), v6xmax(i,j));
                vpymin(i,j) = min(v5ymin(i,j), v6ymin(i,j));
                vpymax(i,j) = max(v5ymax(i,j), v6ymax(i,j));
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Exchange M, Fx, Fy from neighbors
        M = halo_exchange_2d(M, decomp, bc);
        Fx = halo_exchange_2d(Fx, decomp, bc);
        Fy = halo_exchange_2d(Fy, decomp, bc);
        
        % Compute fluxes and wave speeds in halo cells for pas_HLL stencil
        [Fx, Fy, vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext] = ...
            compute_halo_fluxes_and_wavespeeds(M, Fx, Fy, vpxmin, vpxmax, vpymin, vpymax, ...
                                               nx, ny, halo, flag2D_worker, Ma_worker);
        
        % Global reduction for time step (all ranks need same dt)
        vmax_local = max([abs(vpxmax(:)); abs(vpxmin(:)); abs(vpymax(:)); abs(vpymin(:))]);
        vmax = max(gcat(vmax_local, 1));  % Global max across all ranks
        dt = min(CFL_worker*dx_worker/vmax, dtmax_worker);
        dt = min(dt, tmax-t);
        t = t + dt;
        
        % X-direction flux update with processor boundary handling
        Mnpx = apply_flux_update(M, Fx, vpxmin, vpxmax, vpxmin_ext, vpxmax_ext, ...
                                  nx, ny, halo, dt, dx_worker, decomp, 1);
        
        % Y-direction flux update with processor boundary handling
        Mnpy = apply_flux_update(M, Fy, vpymin, vpymax, vpymin_ext, vpymax_ext, ...
                                  nx, ny, halo, dt, dy_worker, decomp, 2);
        
        % Combine updates (Strang splitting) - INTERIOR ONLY
        % Mnpx and Mnpy halos contain stale M values, only interior was updated by pas_HLL
        M(halo+1:halo+nx, halo+1:halo+ny, :) = ...
            Mnpx(halo+1:halo+nx, halo+1:halo+ny, :) + ...
            Mnpy(halo+1:halo+nx, halo+1:halo+ny, :) - ...
            M(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Exchange halos before realizability enforcement
        M = halo_exchange_2d(M, decomp, bc);
        
        % Enforce realizability
        for i = 1:nx
            for j = 1:ny
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                [v6xmin(i,j), v6xmax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D_worker, Ma_worker);
                [v6ymin(i,j), v6ymax(i,j), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D_worker, Ma_worker);
                [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D_worker, Ma_worker);
                Mnp(ih, jh, :) = Mr;
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Apply BGK collision
        for i = 1:nx
            for j = 1:ny
                ih = i + halo;
                jh = j + halo;
                MOM = squeeze(M(ih, jh, :));
                MMC = collision35(MOM, dt, Kn_worker);
                Mnp(ih, jh, :) = MMC;
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :);
        
        % Exchange halos for next iteration
        M = halo_exchange_2d(M, decomp, bc);
        
        % Compute MaxDiff for symmetry check over GLOBAL domain
        % Gather interior data from all ranks to rank 1
        M_interior_local = M(halo+1:halo+nx, halo+1:halo+ny, :);
        
        if spmdIndex == 1
            % Rank 1: collect from all workers to reconstruct global array
            M_global_temp = zeros(Np, Np, Nmom_worker);
            M_global_temp(i0i1(1):i0i1(2), j0j1(1):j0j1(2), :) = M_interior_local;
            
            % Receive from all other ranks
            for src = 2:spmdSize
                data_packet = spmdReceive(src);
                blk = data_packet{1};
                i_range = data_packet{2};
                j_range = data_packet{3};
                M_global_temp(i_range(1):i_range(2), j_range(1):j_range(2), :) = blk;
            end
            
            % Now compute MaxDiff on full global domain
            [~, MaxDiff] = diagnostics('test_symmetry', M_global_temp, Np);
        else
            % All other ranks: send to rank 1
            data_packet = {M_interior_local, i0i1, j0j1};
            spmdSend(data_packet, 1);
            MaxDiff = [];  % Will be broadcast from rank 1
        end
        
        % Compute timing on all ranks
        step_time = toc(step_start_time);
        
        % Compute time per grid point for this rank
        local_grid_points = nx * ny;
        time_per_point_local = step_time / local_grid_points;
        
        % Gather max time per point across all ranks
        max_time_per_point = gop(@max, time_per_point_local);
        
        % Print timestep timing and MaxDiff (only from rank 1)
        if spmdIndex == 1
            fprintf('Step %4d: t = %.6f, dt = %.6e, max s/pt = %.6e s, MaxDiff = %.3e\n', ...
                    nn, t, dt, max_time_per_point, max(abs(MaxDiff)));
        end
    end
    
    % Gather results to rank 1 using asynchronous send/receive
    M_interior = M(halo+1:halo+nx, halo+1:halo+ny, :);
    
    if spmdIndex == 1
        % Rank 1: gather from all ranks using mpi_utils
        M_final = mpi_utils('gather_M', M_interior, i0i1, j0j1, Np, Nmom_worker);
        final_time = t;
        time_steps = nn;
    else
        % All other ranks: send to rank 1 using mpi_utils
        mpi_utils('send_M', M_interior, i0i1, j0j1, 1);
        M_final = [];
        final_time = [];
        time_steps = [];
    end
end

% Extract from composite
if isa(M_final, 'Composite')
    M_final = M_final{1};
    final_time = final_time{1};
    time_steps = time_steps{1};
end

% Grid was created outside spmd, so it's not a Composite
grid_out = grid;

end
