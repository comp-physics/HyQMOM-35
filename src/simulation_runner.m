function [M_final, final_time, time_steps, grid_out] = simulation_runner(params, script_dir)
% Run parallel SPMD simulation loop (3D physical space)
% This is the core time-stepping loop extracted from main.m
%
% Inputs:
%   params     - Struct with all simulation parameters (including Nz)
%   script_dir - Path to project root for worker path setup
%
% Outputs:
%   M_final    - Final moment array (global, gathered to rank 1)
%   final_time - Final simulation time
%   time_steps - Number of time steps taken
%   grid_out   - Grid structure

% MPI parallel execution
spmd
    % Workers need src/ directory in their path
    setup_paths(script_dir);
    
    % Define all constants directly inside spmd block
    halo = 2;  % Required for pas_HLL stencil
    bc = struct('type', 'copy');
    
    % Unpack parameters for workers (avoid Composite overhead in loops)
    Np = params.Np;
    Nz = params.Nz;
    tmax = params.tmax;
    Kn_worker = params.Kn;
    Ma_worker = params.Ma;
    flag2D_worker = params.flag2D;
    CFL_worker = params.CFL;
    dx_worker = params.dx;
    dy_worker = params.dy;
    dz_worker = params.dz;
    Nmom_worker = params.Nmom;
    nnmax_worker = params.nnmax;
    dtmax_worker = params.dtmax;
    
    % Unpack IC parameters (needed for local IC generation)
    rhol_worker = params.rhol;
    rhor_worker = params.rhor;
    T_worker = params.T;
    r110_worker = params.r110;
    r101_worker = params.r101;
    r011_worker = params.r011;
    
    % Unpack diagnostic parameters
    symmetry_check_interval = params.symmetry_check_interval;
    enable_memory_tracking = params.enable_memory_tracking;
    homogeneous_z_worker = params.homogeneous_z;
    
    % Setup domain decomposition (2D in x-y, no decomposition in z)
    decomp = setup_mpi_cartesian_3d(Np, Nz, halo);
    nx = decomp.local_size(1);
    ny = decomp.local_size(2);
    nz = decomp.local_size(3);  % Always equals Nz (no decomposition)
    
    % Create local grid structure
    xmin = -0.5; xmax = 0.5; 
    ymin = -0.5; ymax = 0.5;
    zmin = -0.5; zmax = 0.5;
    dx_global = (xmax - xmin) / Np;
    dy_global = (ymax - ymin) / Np;
    dz_global = (zmax - zmin) / Nz;
    
    % Local grid coordinates (only for this rank's subdomain in x-y)
    i0i1 = decomp.istart_iend;
    j0j1 = decomp.jstart_jend;
    k0k1 = decomp.kstart_kend;  % Always [1, Nz]
    
    % Cell edges for local subdomain
    x_local = xmin + (i0i1(1)-1:i0i1(2)) * dx_global;
    y_local = ymin + (j0j1(1)-1:j0j1(2)) * dy_global;
    z_local = zmin + (k0k1(1)-1:k0k1(2)) * dz_global;
    
    % Cell centers for local subdomain
    xm_local = x_local(1:end-1) + dx_global/2;
    ym_local = y_local(1:end-1) + dy_global/2;
    zm_local = z_local(1:end-1) + dz_global/2;
    
    grid_local = struct('dx', dx_global, 'dy', dy_global, 'dz', dz_global, ...
                        'x', x_local', 'y', y_local', 'z', z_local', ...
                        'xm', xm_local', 'ym', ym_local', 'zm', zm_local');
    
    % Allocate local arrays with halos (4D: x, y, z, moments)
    M = zeros(nx+2*halo, ny+2*halo, nz, Nmom_worker);
    Mnp = M;
    Fx = zeros(nx+2*halo, ny+2*halo, nz, Nmom_worker);
    Fy = zeros(nx+2*halo, ny+2*halo, nz, Nmom_worker);
    Fz = zeros(nx+2*halo, ny+2*halo, nz, Nmom_worker);
    
    % Wave speed arrays (interior only, no halos needed)
    vpxmin = zeros(nx, ny, nz);
    vpxmax = zeros(nx, ny, nz);
    vpymin = zeros(nx, ny, nz);
    vpymax = zeros(nx, ny, nz);
    vpzmin = zeros(nx, ny, nz);
    vpzmax = zeros(nx, ny, nz);
    v5xmin = zeros(nx, ny, nz);
    v5xmax = zeros(nx, ny, nz);
    v5ymin = zeros(nx, ny, nz);
    v5ymax = zeros(nx, ny, nz);
    v5zmin = zeros(nx, ny, nz);
    v5zmax = zeros(nx, ny, nz);
    v6xmin = zeros(nx, ny, nz);
    v6xmax = zeros(nx, ny, nz);
    v6ymin = zeros(nx, ny, nz);
    v6ymax = zeros(nx, ny, nz);
    v6zmin = zeros(nx, ny, nz);
    v6zmax = zeros(nx, ny, nz);
    
    % Build initial conditions locally (no M_global broadcast)
    
    % Background parameters (mean velocities at rest)
    U0 = 0; V0 = 0; W0 = 0;
    
    % Covariance matrix
    C200 = T_worker;
    C020 = T_worker;
    C002 = T_worker;
    C110 = r110_worker * sqrt(C200 * C020);
    C101 = r101_worker * sqrt(C200 * C002);
    C011 = r011_worker * sqrt(C020 * C002);
    
    % Background state (low density)
    Mr_bg = InitializeM4_35(rhor_worker, U0, V0, W0, C200, C110, C101, C020, C011, C002);
    
    % Crossing jets states
    Uc = Ma_worker / sqrt(2);
    Mt = InitializeM4_35(rhol_worker, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002);
    Mb = InitializeM4_35(rhol_worker,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002);
    
    % Jet region bounds (global indices in x-y)
    Csize = floor(0.1 * Np);
    Mint = Np/2 + 1;
    Maxt = Np/2 + 1 + Csize;
    Minb = Np/2 - Csize;
    Maxb = Np/2;
    
    % Fill local subdomain with appropriate IC
    for kk = 1:nz
        gk = k0k1(1) + kk - 1;  % global k index
        z_coord = zm_local(kk);
        
        % Determine if jets exist at this z
        if homogeneous_z_worker
            % Homogeneous in z: jets at all z levels (for validation)
            jets_exist = true;
        else
            % Inhomogeneous: jets only in lower half (z < 0)
            jets_exist = (z_coord < 0.0);
        end
        
        for ii = 1:nx
            gi = i0i1(1) + ii - 1;  % global i index
            for jj = 1:ny
                gj = j0j1(1) + jj - 1;  % global j index
                
                % Default: background
                Mr = Mr_bg;
                
                if jets_exist
                    % Bottom jet (moving up-right)
                    if gi >= Minb && gi <= Maxb && gj >= Minb && gj <= Maxb
                        Mr = Mb;
                    end
                    
                    % Top jet (moving down-left) - overwrites if overlapping
                    if gi >= Mint && gi <= Maxt && gj >= Mint && gj <= Maxt
                        Mr = Mt;
                    end
                end
                
                M(ii + halo, jj + halo, kk, :) = Mr;
            end
        end
    end
    
    % Initial halo exchange
    M = halo_exchange_3d(M, decomp, bc);
    
    % Initialize memory tracking if enabled
    if enable_memory_tracking
        memory_utils('init');
    end
    
    % Time evolution
    t = 0.0;
    nn = 0;
    
    while t < tmax && nn < nnmax_worker
        nn = nn + 1;
        step_start_time = tic;  % Start timing this timestep
        
        % Compute fluxes and wave speeds for interior cells
        for k = 1:nz
            for i = 1:nx
                for j = 1:ny
                    % Access with halo offset
                    ih = i + halo;
                    jh = j + halo;
                    MOM = squeeze(M(ih, jh, k, :));
                    
                    [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                    [v6xmin(i,j,k), v6xmax(i,j,k), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D_worker, Ma_worker);
                    [v6ymin(i,j,k), v6ymax(i,j,k), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D_worker, Ma_worker);
                    [v6zmin(i,j,k), v6zmax(i,j,k), Mr] = eigenvalues6_hyperbolic_3D(Mr, 3, flag2D_worker, Ma_worker);
                    [Mx, My, Mz, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D_worker, Ma_worker);
                    
                    Fx(ih, jh, k, :) = Mx;
                    Fy(ih, jh, k, :) = My;
                    Fz(ih, jh, k, :) = Mz;
                    Mnp(ih, jh, k, :) = Mr;
                    
                    [~, v5xmin(i,j,k), v5xmax(i,j,k)] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
                    [~, v5ymin(i,j,k), v5ymax(i,j,k)] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
                    [~, v5zmin(i,j,k), v5zmax(i,j,k)] = closure_and_eigenvalues(Mr([1,16,20,23,25]));
                    
                    vpxmin(i,j,k) = min(v5xmin(i,j,k), v6xmin(i,j,k));
                    vpxmax(i,j,k) = max(v5xmax(i,j,k), v6xmax(i,j,k));
                    vpymin(i,j,k) = min(v5ymin(i,j,k), v6ymin(i,j,k));
                    vpymax(i,j,k) = max(v5ymax(i,j,k), v6ymax(i,j,k));
                    vpzmin(i,j,k) = min(v5zmin(i,j,k), v6zmin(i,j,k));
                    vpzmax(i,j,k) = max(v5zmax(i,j,k), v6zmax(i,j,k));
                end
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :, :);
        
        % Exchange M, Fx, Fy from neighbors (Fz doesn't need exchange, no z-decomposition)
        M = halo_exchange_3d(M, decomp, bc);
        Fx = halo_exchange_3d(Fx, decomp, bc);
        Fy = halo_exchange_3d(Fy, decomp, bc);
        
        % Compute fluxes and wave speeds in halo cells for pas_HLL stencil
        [Fx, Fy, vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext] = ...
            compute_halo_fluxes_and_wavespeeds_3d(M, Fx, Fy, vpxmin, vpxmax, vpymin, vpymax, vpzmin, vpzmax, ...
                                               nx, ny, nz, halo, flag2D_worker, Ma_worker);
        
        % Global reduction for time step (all ranks need same dt)
        vmax_local = max([abs(vpxmax(:)); abs(vpxmin(:)); abs(vpymax(:)); abs(vpymin(:)); abs(vpzmax(:)); abs(vpzmin(:))]);
        vmax = max(spmdCat(vmax_local, 1));  % Global max across all ranks
        dt = min(CFL_worker*min([dx_worker, dy_worker, dz_worker])/vmax, dtmax_worker);
        dt = min(dt, tmax-t);
        t = t + dt;
        
        % X-direction flux update with processor boundary handling
        Mnpx = apply_flux_update_3d(M, Fx, vpxmin, vpxmax, vpxmin_ext, vpxmax_ext, ...
                                     nx, ny, nz, halo, dt, dx_worker, decomp, 1);
        
        % Y-direction flux update with processor boundary handling
        Mnpy = apply_flux_update_3d(M, Fy, vpymin, vpymax, vpymin_ext, vpymax_ext, ...
                                     nx, ny, nz, halo, dt, dy_worker, decomp, 2);
        
        % Z-direction flux update (simpler: no processor boundaries)
        % Skip if nz=1 (quasi-2D case, no z-transport)
        if nz > 1
            Mnpz = apply_flux_update_3d(M, Fz, vpzmin, vpzmax, [], [], ...
                                         nx, ny, nz, halo, dt, dz_worker, decomp, 3);
        else
            Mnpz = M;  % No z-update for nz=1
        end
        
        % Combine updates (Strang splitting) - INTERIOR ONLY
        % Mnpx, Mnpy, Mnpz halos contain stale M values, only interior was updated by pas_HLL
        if nz > 1
            % Full 3D: x + y + z updates
            M(halo+1:halo+nx, halo+1:halo+ny, :, :) = ...
                Mnpx(halo+1:halo+nx, halo+1:halo+ny, :, :) + ...
                Mnpy(halo+1:halo+nx, halo+1:halo+ny, :, :) + ...
                Mnpz(halo+1:halo+nx, halo+1:halo+ny, :, :) - ...
                2*M(halo+1:halo+nx, halo+1:halo+ny, :, :);
        else
            % Quasi-2D (nz=1): only x + y updates
            M(halo+1:halo+nx, halo+1:halo+ny, :, :) = ...
                Mnpx(halo+1:halo+nx, halo+1:halo+ny, :, :) + ...
                Mnpy(halo+1:halo+nx, halo+1:halo+ny, :, :) - ...
                M(halo+1:halo+nx, halo+1:halo+ny, :, :);
        end
        
        % Exchange halos before realizability enforcement
        M = halo_exchange_3d(M, decomp, bc);
        
        % Enforce realizability
        for k = 1:nz
            for i = 1:nx
                for j = 1:ny
                    ih = i + halo;
                    jh = j + halo;
                    MOM = squeeze(M(ih, jh, k, :));
                    [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D_worker, Ma_worker);
                    [v6xmin(i,j,k), v6xmax(i,j,k), Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D_worker, Ma_worker);
                    [v6ymin(i,j,k), v6ymax(i,j,k), Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D_worker, Ma_worker);
                    [v6zmin(i,j,k), v6zmax(i,j,k), Mr] = eigenvalues6_hyperbolic_3D(Mr, 3, flag2D_worker, Ma_worker);
                    [~,~,~,Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D_worker, Ma_worker);
                    Mnp(ih, jh, k, :) = Mr;
                end
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :, :);
        
        % Apply BGK collision
        for k = 1:nz
            for i = 1:nx
                for j = 1:ny
                    ih = i + halo;
                    jh = j + halo;
                    MOM = squeeze(M(ih, jh, k, :));
                    MMC = collision35(MOM, dt, Kn_worker);
                    Mnp(ih, jh, k, :) = MMC;
                end
            end
        end
        M(halo+1:halo+nx, halo+1:halo+ny, :, :) = Mnp(halo+1:halo+nx, halo+1:halo+ny, :, :);
        
        % Exchange halos for next iteration
        M = halo_exchange_3d(M, decomp, bc);
        
        % Record memory usage for this time step if enabled
        if enable_memory_tracking
            memory_utils('record', nn);
        end
        
        % Compute MaxDiff for symmetry check over GLOBAL domain
        % Only check every symmetry_check_interval steps to reduce overhead
        % OPTIMIZATION: Gather only diagonal entries (not full field) to minimize memory
        if mod(nn, symmetry_check_interval) == 0 || nn == 1
            % Find global diagonal indices present on this rank
            % Diagonal: global indices where i == j, use middle z-plane
            kmid = floor(nz/2) + 1;
            diag_i_global = max(i0i1(1), j0j1(1)) : min(i0i1(2), j0j1(2));
            ndiag = numel(diag_i_global);
            
            % Extract the 5 moment components needed for symmetry test from diagonal cells
            diag5_local = zeros(ndiag, 5);
            for idx = 1:ndiag
                gi = diag_i_global(idx);
                ii = gi - i0i1(1) + 1;  % local i index (interior, 1-based)
                jj = gi - j0j1(1) + 1;  % local j index (interior, 1-based)
                ih = ii + halo;
                jh = jj + halo;
                diag5_local(idx, 1) = M(ih, jh, kmid, 1);
                diag5_local(idx, 2) = M(ih, jh, kmid, 2);
                diag5_local(idx, 3) = M(ih, jh, kmid, 3);
                diag5_local(idx, 4) = M(ih, jh, kmid, 4);
                diag5_local(idx, 5) = M(ih, jh, kmid, 5);
            end
            
            if spmdIndex == 1
                % Rank 1: assemble full diagonal from all ranks
                diag5_global = zeros(Np, 5);
                diag5_global(diag_i_global, :) = diag5_local;
                
                % Receive diagonal slices from other ranks
                for src = 2:spmdSize
                    pkt = spmdReceive(src);  % {gi_range, diag5_vals}
                    gi_range = pkt{1};
                    vals = pkt{2};
                    diag5_global(gi_range, :) = vals;
                end
                
                % Compute Diff and MaxDiff using only diagonal entries
                % This is the core of test_symmetry: M(i,i,k) vs M(Np+1-i,Np+1-i,k)
                Diff = zeros(Np, 5);
                for i = 1:Np
                    j = Np + 1 - i;
                    Diff(i, 1) = diag5_global(i, 1) - diag5_global(j, 1);
                    Diff(i, 2) = diag5_global(i, 2) + diag5_global(j, 2);
                    Diff(i, 3) = diag5_global(i, 3) - diag5_global(j, 3);
                    Diff(i, 4) = diag5_global(i, 4) + diag5_global(j, 4);
                    Diff(i, 5) = diag5_global(i, 5) - diag5_global(j, 5);
                end
                MaxDiff = zeros(5, 1);
                for k = 1:5
                    Normk = norm(Diff(:, k));
                    MaxDiff(k) = max(Diff(:, k)) / (Normk + 1);
                end
            else
                % Non-root ranks: send diagonal slice to rank 1
                spmdSend({diag_i_global, diag5_local}, 1);
                MaxDiff = [];
            end
        else
            % Skip symmetry check this step
            MaxDiff = 0;
        end
        
        % Compute timing on all ranks
        step_time = toc(step_start_time);
        
        % Compute time per grid point for this rank
        local_grid_points = nx * ny * nz;
        time_per_point_local = step_time / local_grid_points;
        
        % Gather max time per point across all ranks
        max_time_per_point = spmdReduce(@max, time_per_point_local);
        
        % Print timestep timing and MaxDiff (only from rank 1)
        if spmdIndex == 1
            if mod(nn, symmetry_check_interval) == 0 || nn == 1
                fprintf('Step %4d: t = %.6f, dt = %.6e, max s/pt = %.6e s, MaxDiff = %.3e\n', ...
                        nn, t, dt, max_time_per_point, max(abs(MaxDiff)));
            else
                fprintf('Step %4d: t = %.6f, dt = %.6e, max s/pt = %.6e s\n', ...
                        nn, t, dt, max_time_per_point);
            end
        end
    end
    
    % Report memory usage statistics if enabled
    if enable_memory_tracking
        memory_utils('report', spmdSize);
    end
    
    % Gather results to rank 1 using asynchronous send/receive
    M_interior = M(halo+1:halo+nx, halo+1:halo+ny, :, :);
    
    if spmdIndex == 1
        % Rank 1: gather from all ranks using mpi_utils
        M_final = mpi_utils('gather_M_3d', M_interior, i0i1, j0j1, k0k1, Np, Nz, Nmom_worker);
        final_time = t;
        time_steps = nn;
        
        % Rank 1 creates the global grid structure for output
        grid_out_local = struct();
        grid_out_local.x = linspace(xmin, xmax, Np+1)';
        grid_out_local.y = linspace(ymin, ymax, Np+1)';
        grid_out_local.z = linspace(zmin, zmax, Nz+1)';
        grid_out_local.dx = dx_global;
        grid_out_local.dy = dy_global;
        grid_out_local.dz = dz_global;
        grid_out_local.xm = grid_out_local.x(1:Np) + dx_global/2;
        grid_out_local.ym = grid_out_local.y(1:Np) + dy_global/2;
        grid_out_local.zm = grid_out_local.z(1:Nz) + dz_global/2;
    else
        % All other ranks: send to rank 1 using mpi_utils
        mpi_utils('send_M_3d', M_interior, i0i1, j0j1, k0k1, 1);
        M_final = [];
        final_time = [];
        time_steps = [];
        grid_out_local = [];
    end
end

% Extract from composite
if isa(M_final, 'Composite')
    M_final = M_final{1};
    final_time = final_time{1};
    time_steps = time_steps{1};
    grid_out = grid_out_local{1};
else
    grid_out = grid_out_local;
end

end
