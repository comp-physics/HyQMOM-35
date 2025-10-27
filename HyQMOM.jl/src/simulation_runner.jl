using LinearAlgebra: norm
using JLD2
using Printf

"""
    simulation_runner(params)

Run MPI-parallel simulation loop.

This is the core time-stepping loop that orchestrates all components:
- Domain decomposition
- Initial conditions
- Time evolution with flux computation
- Halo exchange
- Realizability enforcement
- BGK collision
- Result gathering

# Arguments
- `params`: Named tuple with simulation parameters:
  - `Nx`: Global grid size in x direction
  - `Ny`: Global grid size in y direction
  - `Nz`: Global grid size in z direction
  - `tmax`: Maximum simulation time
  - `Kn`: Knudsen number
  - `Ma`: Mach number
  - `flag2D`: 2D flag (1 for 2D, 0 for 3D)
  - `CFL`: CFL number
  - `dx`, `dy`, `dz`: Grid spacing
  - `Nmom`: Number of moments (35)
  - `nnmax`: Maximum number of time steps
  - `dtmax`: Maximum time step size
  - IC parameters: `rhol`, `rhor`, `T`, `r110`, `r101`, `r011`
  - `symmetry_check_interval`: How often to check symmetry
  - `homogeneous_z`: Whether jets exist at all z levels (validation mode)
  - `enable_memory_tracking`: Enable memory tracking (not implemented)
  - `snapshot_interval`: (optional) Save snapshots every N steps. If not provided or 0, no snapshots saved.
  - Standardized moments (S4) and central moments (C4) are automatically saved with each snapshot

# Returns
If `snapshot_interval` is provided and > 0:
- Snapshots are streamed directly to a JLD2 file (one snapshot at a time, never all in memory)
- On rank 0: `snapshot_filename` (path to JLD2 file), `grid_out`
- On other ranks: `nothing`, `nothing`
  
Snapshot file structure:
- `meta/params`: Simulation parameters
- `meta/snapshot_interval`: Interval value
- `meta/n_snapshots`: Total number of snapshots saved
- `grid`: Grid structure
- `snapshots/000001/{M, t, step[, S][, C]}`: First snapshot
- `snapshots/000002/...`: Second snapshot, etc.

Otherwise (standard mode):
- `M_final`: Final moment array (global, gathered to rank 0), or `nothing` on other ranks
- `final_time`: Final simulation time
- `time_steps`: Number of time steps taken
- `grid_out`: Grid structure (only on rank 0)

# Algorithm
See simulation_runner.m for detailed algorithm description.
"""
function simulation_runner(params)
    # MPI initialization should be done before calling this function
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    # Constants
    halo = 2  # Required for pas_HLL stencil
    bc = :copy
    
    # Unpack parameters
    Nx = params.Nx
    Ny = params.Ny
    Nz = params.Nz
    tmax = params.tmax
    Kn = params.Kn
    Ma = params.Ma
    flag2D = params.flag2D
    CFL = params.CFL
    Nmom = params.Nmom
    nnmax = params.nnmax
    dtmax = params.dtmax
    
    # IC parameters
    rhol = params.rhol
    rhor = params.rhor
    T = params.T
    r110 = params.r110
    r101 = params.r101
    r011 = params.r011
    
    # Diagnostic parameters
    symmetry_check_interval = params.symmetry_check_interval
    homogeneous_z = params.homogeneous_z
    debug_output = params.debug_output
    
    # Snapshot saving parameters
    snapshot_interval = get(params, :snapshot_interval, 0)
    save_snapshots = (snapshot_interval > 0)
    save_standardized = true  # Always save S4 with snapshots (required for moment space visualization)
    save_central = true  # Always save C4 with snapshots
    
    # Streaming snapshot file (rank 0 only)
    snapshot_filename = get(params, :snapshot_filename, nothing)
    jld_file = nothing
    snap_count = 0
    
    # Setup domain decomposition (2D in x-y, no decomposition in z)
    decomp = setup_mpi_cartesian_3d(Nx, Ny, Nz, halo, comm)
    nx = decomp.local_size[1]
    ny = decomp.local_size[2]
    nz = decomp.local_size[3]  # Always equals Nz (no decomposition)
    
    # Domain extents (allow user to specify, default to [-0.5, 0.5])
    xmin = get(params, :xmin, -0.5)
    xmax = get(params, :xmax,  0.5)
    ymin = get(params, :ymin, -0.5)
    ymax = get(params, :ymax,  0.5)
    zmin = get(params, :zmin, -0.5)
    zmax = get(params, :zmax,  0.5)
    
    # Compute grid spacing from domain extents and resolution
    # Note: dx, dy, dz are always computed automatically and never read from params
    dx_global = (xmax - xmin) / Nx
    dy_global = (ymax - ymin) / Ny
    dz_global = (zmax - zmin) / Nz
    
    dx = dx_global
    dy = dy_global
    dz = dz_global
    
    # Local grid coordinates (only for this rank's subdomain in x-y)
    i0i1 = decomp.istart_iend
    j0j1 = decomp.jstart_jend
    k0k1 = decomp.kstart_kend  # Always (1, Nz)
    
    # Cell edges for local subdomain
    x_local = range(xmin + (i0i1[1]-1)*dx_global, step=dx_global, length=nx+1)
    y_local = range(ymin + (j0j1[1]-1)*dy_global, step=dy_global, length=ny+1)
    z_local = range(zmin + (k0k1[1]-1)*dz_global, step=dz_global, length=nz+1)
    
    # Cell centers
    xm_local = collect(x_local[1:end-1]) .+ dx_global/2
    ym_local = collect(y_local[1:end-1]) .+ dy_global/2
    zm_local = collect(z_local[1:end-1]) .+ dz_global/2
    
    grid_local = (dx=dx_global, dy=dy_global, dz=dz_global,
                  x=collect(x_local), y=collect(y_local), z=collect(z_local),
                  xm=xm_local, ym=ym_local, zm=zm_local)
    
    # Allocate local arrays with halos (4D: x, y, z, moments)
    M = zeros(Float64, nx+2*halo, ny+2*halo, nz, Nmom)
    Mnp = similar(M)
    Fx = zeros(Float64, nx+2*halo, ny+2*halo, nz, Nmom)
    Fy = zeros(Float64, nx+2*halo, ny+2*halo, nz, Nmom)
    Fz = zeros(Float64, nx+2*halo, ny+2*halo, nz, Nmom)
    
    # Wave speed arrays (interior only, no halos needed)
    vpxmin = zeros(Float64, nx, ny, nz)
    vpxmax = zeros(Float64, nx, ny, nz)
    vpymin = zeros(Float64, nx, ny, nz)
    vpymax = zeros(Float64, nx, ny, nz)
    vpzmin = zeros(Float64, nx, ny, nz)
    vpzmax = zeros(Float64, nx, ny, nz)
    v5xmin = zeros(Float64, nx, ny, nz)
    v5xmax = zeros(Float64, nx, ny, nz)
    v5ymin = zeros(Float64, nx, ny, nz)
    v5ymax = zeros(Float64, nx, ny, nz)
    v5zmin = zeros(Float64, nx, ny, nz)
    v5zmax = zeros(Float64, nx, ny, nz)
    v6xmin = zeros(Float64, nx, ny, nz)
    v6xmax = zeros(Float64, nx, ny, nz)
    v6ymin = zeros(Float64, nx, ny, nz)
    v6ymax = zeros(Float64, nx, ny, nz)
    v6zmin = zeros(Float64, nx, ny, nz)
    v6zmax = zeros(Float64, nx, ny, nz)
    
    # Build initial conditions locally
    # Check if custom ICs are provided
    if haskey(params, :use_custom_ic) && params.use_custom_ic
        # Use flexible custom IC system
        # Need global grid for custom ICs
        xm_global = collect(range(xmin + dx_global/2, step=dx_global, length=Nx))
        ym_global = collect(range(ymin + dy_global/2, step=dy_global, length=Ny))
        zm_global = collect(range(zmin + dz_global/2, step=dz_global, length=Nz))
        
        grid_params = (Nx=Nx, Ny=Ny, Nz=Nz, xm=xm_global, ym=ym_global, zm=zm_global)
        M = initialize_moment_field_mpi(decomp, grid_params, 
                                       params.ic_background, params.ic_jets;
                                       r110=r110, r101=r101, r011=r011, halo=halo)
    else
        # Use original hardcoded crossing jets IC
        U0, V0, W0 = 0.0, 0.0, 0.0
        
        # Covariance matrix
        C200 = T
        C020 = T
        C002 = T
        C110 = r110 * sqrt(C200 * C020)
        C101 = r101 * sqrt(C200 * C002)
        C011 = r011 * sqrt(C020 * C002)
        
        # Background state (low density)
        Mr_bg = InitializeM4_35(rhor, U0, V0, W0, C200, C110, C101, C020, C011, C002)
        
        # Crossing jets states
        Uc = Ma / sqrt(2.0)
        Mt = InitializeM4_35(rhol, -Uc, -Uc, W0, C200, C110, C101, C020, C011, C002)
        Mb = InitializeM4_35(rhol,  Uc,  Uc, W0, C200, C110, C101, C020, C011, C002)
        
        # Jet region bounds (global indices)
        # x-y plane: 10% of domain size
        Csize_x = floor(Int, 0.1 * Nx)
        Csize_y = floor(Int, 0.1 * Ny)
        Mint_x = div(Nx, 2) + 1
        Maxt_x = div(Nx, 2) + 1 + Csize_x
        Minb_x = div(Nx, 2) - Csize_x
        Maxb_x = div(Nx, 2)
        Mint_y = div(Ny, 2) + 1
        Maxt_y = div(Ny, 2) + 1 + Csize_y
        Minb_y = div(Ny, 2) - Csize_y
        Maxb_y = div(Ny, 2)
        
        # z-direction: Create cubes (not extruded squares)
        # Cubes centered at z=0
        Csize_z = min(Csize_x, Csize_y)  # Use smaller of x,y sizes for cubic regions
        Mint_z = div(Nz, 2) + 1 - div(Csize_z, 2)
        Maxt_z = div(Nz, 2) + 1 + div(Csize_z, 2)
        
        # Fill local subdomain with appropriate IC (sharp transitions)
        # Note: Sharp transitions work better than smooth blending because
        # blending moment vectors doesn't preserve realizability
        for kk in 1:nz
            gk = k0k1[1] + kk - 1  # global k index
            
            for ii in 1:nx
                gi = i0i1[1] + ii - 1  # global i index
                for jj in 1:ny
                    gj = j0j1[1] + jj - 1  # global j index
                    
                    # Default: background
                    Mr = Mr_bg
                    
                    # Bottom jet cube: check x, y, AND z bounds
                    if (gi >= Minb_x && gi <= Maxb_x && 
                        gj >= Minb_y && gj <= Maxb_y &&
                        gk >= Mint_z && gk <= Maxt_z)
                        Mr = Mb
                    end
                    
                    # Top jet cube: check x, y, AND z bounds (overwrites if overlapping)
                    if (gi >= Mint_x && gi <= Maxt_x && 
                        gj >= Mint_y && gj <= Maxt_y &&
                        gk >= Mint_z && gk <= Maxt_z)
                        Mr = Mt
                    end
                    
                    M[ii + halo, jj + halo, kk, :] = Mr
                end
            end
        end
    end
    
    # Initial halo exchange
    halo_exchange_3d!(M, decomp, bc)
    
    # Save initial condition as first snapshot if requested
    i0i1 = decomp.istart_iend
    j0j1 = decomp.jstart_jend
    k0k1 = decomp.kstart_kend
    
    # Initialize streaming snapshot file (rank 0 only)
    if save_snapshots && rank == 0
        if snapshot_filename === nothing
            # Auto-generate filename with simulation parameters
            snapshot_filename = @sprintf(
                "snapshots_Nx%d_Ny%d_Nz%d_Kn%.2f_Ma%.2f_tmax%.4f.jld2",
                Nx, Ny, Nz, Kn, Ma, tmax
            )
        end
        
        jld_file = jldopen(snapshot_filename, "w")
        jld_file["meta/params"] = params
        jld_file["meta/snapshot_interval"] = snapshot_interval
        # Don't write n_snapshots yet - will write final count at end
        
        println("  Streaming snapshots to: $snapshot_filename")
    end
    
    if save_snapshots
        M_interior = M[halo+1:halo+nx, halo+1:halo+ny, :, :]
        if rank == 0
            M_gathered = gather_M(M_interior, i0i1, j0j1, k0k1, Nx, Ny, Nz, Nmom, comm)
            
            # Write snapshot to file (streaming mode)
            snap_count += 1
            snap_key = lpad(snap_count, 6, '0')
            jld_file["snapshots/$snap_key/M"] = M_gathered
            jld_file["snapshots/$snap_key/t"] = 0.0
            jld_file["snapshots/$snap_key/step"] = 0
            
            if save_standardized
                S = compute_standardized_field(M_gathered)
                jld_file["snapshots/$snap_key/S"] = S
            end
            if save_central
                C = compute_central_field(M_gathered)
                jld_file["snapshots/$snap_key/C"] = C
            end
            
            # Free memory immediately
            M_gathered = nothing
            GC.gc()
            
            if debug_output
                println("Streamed snapshot 1: t=0.0, step=0")
            end
        else
            gather_M(M_interior, i0i1, j0j1, k0k1, Nx, Ny, Nz, Nmom, comm)
        end
    end
    
    # Time evolution
    t = 0.0
    nn = 0
    
    if rank == 0
        println("Starting time evolution...")
        println("  Grid: $(Nx)x$(Ny)x$(Nz), Ranks: $(nprocs), Local: $(nx)x$(ny)x$(nz)")
        println("  tmax: $(tmax), CFL: $(CFL), Ma: $(Ma), Kn: $(Kn)")
        println("  homogeneous_z: $(homogeneous_z)")
        if save_snapshots
            println("  Snapshot interval: every $(snapshot_interval) steps")
        end
    end
    
    while t < tmax && nn < nnmax
        nn += 1
        step_start_time = time()
        
        
        # Compute fluxes and wave speeds for interior cells
        for k in 1:nz
            for i in 1:nx
                for j in 1:ny
                    ih = i + halo
                    jh = j + halo
                    MOM = M[ih, jh, k, :]
                    
                    _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                    v6xmin[i,j,k], v6xmax[i,j,k], Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                    v6ymin[i,j,k], v6ymax[i,j,k], Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                    v6zmin[i,j,k], v6zmax[i,j,k], Mr = eigenvalues6z_hyperbolic_3D(Mr, flag2D, Ma)
                    Mx, My, Mz, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
                    
                    Fx[ih, jh, k, :] = Mx
                    Fy[ih, jh, k, :] = My
                    Fz[ih, jh, k, :] = Mz
                    Mnp[ih, jh, k, :] = Mr
                    
                    _, v5xmin[i,j,k], v5xmax[i,j,k] = closure_and_eigenvalues(Mr[[1,2,3,4,5]])
                    _, v5ymin[i,j,k], v5ymax[i,j,k] = closure_and_eigenvalues(Mr[[1,6,10,13,15]])
                    _, v5zmin[i,j,k], v5zmax[i,j,k] = closure_and_eigenvalues(Mr[[1,16,20,23,25]])
                    
                    vpxmin[i,j,k] = min(v5xmin[i,j,k], v6xmin[i,j,k])
                    vpxmax[i,j,k] = max(v5xmax[i,j,k], v6xmax[i,j,k])
                    vpymin[i,j,k] = min(v5ymin[i,j,k], v6ymin[i,j,k])
                    vpymax[i,j,k] = max(v5ymax[i,j,k], v6ymax[i,j,k])
                    vpzmin[i,j,k] = min(v5zmin[i,j,k], v6zmin[i,j,k])
                    vpzmax[i,j,k] = max(v5zmax[i,j,k], v6zmax[i,j,k])
                end
            end
        end
        M[halo+1:halo+nx, halo+1:halo+ny, :, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :, :]
        
        # Exchange M, Fx, Fy, Fz from neighbors
        halo_exchange_3d!(M, decomp, bc)
        halo_exchange_3d!(Fx, decomp, bc)
        halo_exchange_3d!(Fy, decomp, bc)
        halo_exchange_3d!(Fz, decomp, bc)
        
        # Compute fluxes and wave speeds in halo cells
        vpxmin_ext = zeros(Float64, nx+2*halo, ny, nz)
        vpxmax_ext = zeros(Float64, nx+2*halo, ny, nz)
        vpymin_ext = zeros(Float64, nx, ny+2*halo, nz)
        vpymax_ext = zeros(Float64, nx, ny+2*halo, nz)
        
        compute_halo_fluxes_and_wavespeeds_3d!(Fx, Fy, vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext,
                                               M, vpxmin, vpxmax, vpymin, vpymax, vpzmin, vpzmax,
                                               nx, ny, nz, halo, flag2D, Ma)
        
        # Global reduction for time step
        vmax_local = maximum([abs.(vpxmax); abs.(vpxmin); abs.(vpymax); abs.(vpymin); abs.(vpzmax); abs.(vpzmin)])
        vmax = MPI.Allreduce(vmax_local, max, comm)
        dt = min(CFL*min(dx,dy,dz)/vmax, dtmax)
        dt = min(dt, tmax-t)
        
        # Debug: print vmax if it's unusually large
        if rank == 0 && (vmax > 100.0 || isnan(vmax) || isinf(vmax))
            @printf("  WARNING: vmax = %.6e (unusually large or invalid)\n", vmax)
            # Find which cell has the problem (only check first z-slice for brevity)
            for i in 1:nx
                for j in 1:ny
                    k = 1
                    if abs(vpxmax[i,j,k]) > 100.0 || abs(vpxmin[i,j,k]) > 100.0 ||
                       abs(vpymax[i,j,k]) > 100.0 || abs(vpymin[i,j,k]) > 100.0 ||
                       abs(vpzmax[i,j,k]) > 100.0 || abs(vpzmin[i,j,k]) > 100.0 ||
                       isnan(vpxmax[i,j,k]) || isnan(vpxmin[i,j,k]) ||
                       isnan(vpymax[i,j,k]) || isnan(vpymin[i,j,k]) ||
                       isnan(vpzmax[i,j,k]) || isnan(vpzmin[i,j,k])
                        @printf("    Cell (%d,%d,%d): vpxmin=%.2e, vpxmax=%.2e, vpymin=%.2e, vpymax=%.2e, vpzmin=%.2e, vpzmax=%.2e\n",
                               i, j, k, vpxmin[i,j,k], vpxmax[i,j,k], vpymin[i,j,k], vpymax[i,j,k], vpzmin[i,j,k], vpzmax[i,j,k])
                        # Print moments at this cell
                        ih = i + halo
                        jh = j + halo
                        @printf("    M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n",
                               M[ih,jh,k,1], M[ih,jh,k,2], M[ih,jh,k,3], M[ih,jh,k,4], M[ih,jh,k,5])
                    end
                end
            end
        end
        
        t += dt
        
        # X-direction flux update
        Mnpx = similar(M)
        apply_flux_update_3d!(Mnpx, M, Fx, vpxmin, vpxmax, vpxmin_ext, vpxmax_ext,
                              nx, ny, nz, halo, dt, dx, decomp, 1)
        
        # Y-direction flux update
        Mnpy = similar(M)
        apply_flux_update_3d!(Mnpy, M, Fy, vpymin, vpymax, vpymin_ext, vpymax_ext,
                              nx, ny, nz, halo, dt, dy, decomp, 2)
        
        # Z-direction flux update
        Mnpz = similar(M)
        vpzmin_ext = zeros(Float64, nx+2*halo, ny, nz)  # Not used for Z (no halo extension)
        vpzmax_ext = zeros(Float64, nx+2*halo, ny, nz)  # Not used for Z (no halo extension)
        apply_flux_update_3d!(Mnpz, M, Fz, vpzmin, vpzmax, vpzmin_ext, vpzmax_ext,
                              nx, ny, nz, halo, dt, dz, decomp, 3)
        
        # Combine updates (Strang splitting) - INTERIOR ONLY
        M[halo+1:halo+nx, halo+1:halo+ny, :, :] =
            Mnpx[halo+1:halo+nx, halo+1:halo+ny, :, :] +
            Mnpy[halo+1:halo+nx, halo+1:halo+ny, :, :] +
            Mnpz[halo+1:halo+nx, halo+1:halo+ny, :, :] -
            2.0 .* M[halo+1:halo+nx, halo+1:halo+ny, :, :]
        
        # Exchange halos before realizability enforcement
        halo_exchange_3d!(M, decomp, bc)
        
        # Enforce realizability
        for k in 1:nz
            for i in 1:nx
                for j in 1:ny
                    ih = i + halo
                    jh = j + halo
                    MOM = M[ih, jh, k, :]
                    
                    _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                    v6xmin[i,j,k], v6xmax[i,j,k], Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                    v6ymin[i,j,k], v6ymax[i,j,k], Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                    _, _, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
                    
                    Mnp[ih, jh, k, :] = Mr
                end
            end
        end
        
        M[halo+1:halo+nx, halo+1:halo+ny, :, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :, :]
        
        # Apply BGK collision
        for k in 1:nz
            for i in 1:nx
                for j in 1:ny
                    ih = i + halo
                    jh = j + halo
                    MOM = M[ih, jh, k, :]
                    MMC = collision35(MOM, dt, Kn)
                    Mnp[ih, jh, k, :] = MMC
                end
            end
        end
        M[halo+1:halo+nx, halo+1:halo+ny, :, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :, :]
        
        # Exchange halos for next iteration
        halo_exchange_3d!(M, decomp, bc)
        
        # Symmetry check: test that M(i,i,k) ≈ M(Nx+1-i, Ny+1-i, k) along diagonal (for square grids)
        # Only check every symmetry_check_interval steps to reduce overhead
        MaxDiff = 0.0
        if mod(nn, symmetry_check_interval) == 0 || nn == 1
            # Gather diagonal entries from all ranks
            # Diagonal: global indices where i == j
            i0, i1 = decomp.istart_iend
            j0, j1 = decomp.jstart_jend
            i_start = max(i0, j0)
            i_end = min(i1, j1)
            diag_i_global = i_start:i_end
            ndiag = length(diag_i_global)
            
            # Extract 5 moment components from diagonal cells for symmetry test
            # MATLAB checks indices 1-5: M000, M100, M200, M300, M400
            diag5_local = zeros(ndiag, 5)
            for (idx, gi) in enumerate(diag_i_global)
                ii = gi - i0 + 1  # local i index (interior, 1-based)
                jj = gi - j0 + 1  # local j index (interior, 1-based)
                diag5_local[idx, 1] = M[halo+ii, halo+jj, 1, 1]   # M000 (first z-slice)
                diag5_local[idx, 2] = M[halo+ii, halo+jj, 1, 2]   # M100
                diag5_local[idx, 3] = M[halo+ii, halo+jj, 1, 3]   # M200
                diag5_local[idx, 4] = M[halo+ii, halo+jj, 1, 4]   # M300
                diag5_local[idx, 5] = M[halo+ii, halo+jj, 1, 5]   # M400
            end
            
            # Gather to rank 0 (only works for square grids where Nx==Ny)
            if rank == 0
                diag5_global = zeros(min(Nx, Ny), 5)
                # Copy local contribution
                for (idx, gi) in enumerate(diag_i_global)
                    diag5_global[gi, :] = diag5_local[idx, :]
                end
                # Receive from other ranks
                for src in 1:(nprocs-1)
                    ndiag_remote_buf = Vector{Int}(undef, 1)
                    MPI.Recv!(ndiag_remote_buf, comm; source=src, tag=0)
                    ndiag_remote = ndiag_remote_buf[1]
                    if ndiag_remote > 0
                        i_remote = Vector{Int}(undef, ndiag_remote)
                        MPI.Recv!(i_remote, comm; source=src, tag=1)
                        diag_remote = Array{Float64}(undef, ndiag_remote, 5)
                        MPI.Recv!(diag_remote, comm; source=src, tag=2)
                        for idx in 1:ndiag_remote
                            diag5_global[i_remote[idx], :] = diag_remote[idx, :]
                        end
                    end
                end
                
                # Compute symmetry differences: M(i,i,k) vs M(Nx+1-i,Ny+1-i,k)
                Ndiag = min(Nx, Ny)
                Diff = zeros(Ndiag, 5)
                for i = 1:Ndiag
                    j = Ndiag + 1 - i
                    Diff[i, 1] = diag5_global[i, 1] - diag5_global[j, 1]  # M000 (symmetric)
                    Diff[i, 2] = diag5_global[i, 2] + diag5_global[j, 2]  # M100 (antisymmetric)
                    Diff[i, 3] = diag5_global[i, 3] - diag5_global[j, 3]  # M010 (symmetric)
                    Diff[i, 4] = diag5_global[i, 4] + diag5_global[j, 4]  # M001 (antisymmetric)
                    Diff[i, 5] = diag5_global[i, 5] - diag5_global[j, 5]  # M002 (symmetric)
                end
                # Compute normalized MaxDiff for each moment (matching MATLAB)
                MaxDiff_vec = zeros(5)
                for k = 1:5
                    Normk = norm(Diff[:, k])
                    MaxDiff_vec[k] = maximum(Diff[:, k]) / (Normk + 1)
                end
                MaxDiff = maximum(abs.(MaxDiff_vec))
            else
                # Send local diagonal to rank 0
                MPI.Send([ndiag], comm; dest=0, tag=0)
                if ndiag > 0
                    MPI.Send(collect(diag_i_global), comm; dest=0, tag=1)
                    MPI.Send(diag5_local, comm; dest=0, tag=2)
                end
            end
        end
        
        # Timing: compute per global grid point for fair comparison across ranks
        step_time = time() - step_start_time
        max_step_time = MPI.Allreduce(step_time, max, comm)
        global_grid_points = Nx * Ny * Nz
        time_per_point_global = max_step_time / global_grid_points
        
        # Print timestep info
        if rank == 0
            if mod(nn, symmetry_check_interval) == 0 || nn == 1
                @printf("Step %4d: t = %.4f, dt = %.4e, wall = %.4f s, s/pt = %.4e s, MaxDiff = %.3e\n",
                       nn, t, dt, max_step_time, time_per_point_global, MaxDiff)
            else
                @printf("Step %4d: t = %.4f, dt = %.4e, wall = %.4f s, s/pt = %.4e s\n",
                       nn, t, dt, max_step_time, time_per_point_global)
            end
        end
        
        # Save snapshot if requested and it's time
        if save_snapshots && (mod(nn, snapshot_interval) == 0 || t >= tmax || nn == nnmax)
            M_interior = M[halo+1:halo+nx, halo+1:halo+ny, :, :]
            if rank == 0
                M_gathered = gather_M(M_interior, i0i1, j0j1, k0k1, Nx, Ny, Nz, Nmom, comm)
                
                # Write snapshot to file (streaming mode)
                snap_count += 1
                snap_key = lpad(snap_count, 6, '0')
                jld_file["snapshots/$snap_key/M"] = M_gathered
                jld_file["snapshots/$snap_key/t"] = t
                jld_file["snapshots/$snap_key/step"] = nn
                
                if save_standardized
                    S = compute_standardized_field(M_gathered)
                    jld_file["snapshots/$snap_key/S"] = S
                end
                if save_central
                    C = compute_central_field(M_gathered)
                    jld_file["snapshots/$snap_key/C"] = C
                end
                
                # Free memory immediately
                M_gathered = nothing
                GC.gc()
                
                println("Streamed snapshot $snap_count: t=$(round(t, digits=4)), step=$nn")
            else
                gather_M(M_interior, i0i1, j0j1, k0k1, Nx, Ny, Nz, Nmom, comm)
            end
        end
    end
    
    if rank == 0
        println("Time evolution complete: $(nn) steps, t = $(t)")
        if save_snapshots
            println("Streamed $snap_count snapshots total")
        end
    end
    
    # Create global grid structure
    grid_out = nothing
    if rank == 0
        grid_out = (x = collect(range(xmin, xmax, length=Nx+1)),
                   y = collect(range(ymin, ymax, length=Ny+1)),
                   z = collect(range(zmin, zmax, length=Nz+1)),
                   dx = dx_global,
                   dy = dy_global,
                   dz = dz_global,
                   xm = collect(range(xmin + dx_global/2, step=dx_global, length=Nx)),
                   ym = collect(range(ymin + dy_global/2, step=dy_global, length=Ny)),
                   zm = collect(range(zmin + dz_global/2, step=dz_global, length=Nz)))
    end
    
    # Return depends on whether snapshots were requested
    if save_snapshots
        # Streaming snapshot mode: close file and return filename
        if rank == 0
            jld_file["grid"] = grid_out
            jld_file["meta/n_snapshots"] = snap_count  # Write final count
            close(jld_file)
            println("Snapshot file closed: $snapshot_filename")
            return snapshot_filename, grid_out
        else
            return nothing, nothing
        end
    else
        # Standard mode: gather final result
        M_interior = M[halo+1:halo+nx, halo+1:halo+ny, :, :]
        
        if rank == 0
            M_final = gather_M(M_interior, i0i1, j0j1, k0k1, Nx, Ny, Nz, Nmom, comm)
            final_time = t
            time_steps = nn
            return M_final, final_time, time_steps, grid_out
        else
            # Other ranks send to rank 0
            gather_M(M_interior, i0i1, j0j1, k0k1, Nx, Ny, Nz, Nmom, comm)
            return nothing, t, nn, nothing
        end
    end
end

"""
    run_simulation(; Np=nothing, Nx=20, Ny=20, tmax=0.1, num_workers=1, kwargs...)

High-level wrapper for running HyQMOM simulations with sensible defaults.

This function provides a convenient interface to the HyQMOM solver with automatic
MPI initialization, parameter setup, and optional visualization. It handles both
serial and parallel execution transparently.

# Arguments
- `Np`: (Legacy) Global grid size for square grid (Np x Np), overrides Nx/Ny if provided
- `Nx`: Global grid size in x direction, default 20
- `Ny`: Global grid size in y direction, default 20  
- `Nz`: Global grid size in z direction, default 1
- `tmax`: Maximum simulation time, default 0.1
- `num_workers`: Number of MPI ranks (must match mpiexec -n), default 1
- `verbose`: Print progress information, default true
- `save_output`: Save output to file, default false
- `enable_plots`: Generate PyPlot visualization after simulation, default false
- `save_figures`: Save figures to disk (requires enable_plots=true), default false
- `output_dir`: Directory for saved figures, default "."
- `Kn`: Knudsen number (controls collision frequency), default 1.0
- `Ma`: Mach number (controls jet velocity), default 0.0
- `flag2D`: 2D flag (1 for 2D, 0 for 3D), default 0
- `CFL`: CFL number for numerical stability, default 0.5
- `homogeneous_z`: Jets at all z levels (true) or only lower half (false), default true

# Returns
Dictionary with simulation results:
- `:M`: Final moment field (Nx x Ny x Nz x 35) on rank 0, nothing on other ranks
- `:final_time`: Actual final time reached
- `:time_steps`: Number of time steps taken
- `:xm`, `:ym`, `:zm`: Grid coordinates (only on rank 0)
- `:Nx`, `:Ny`, `:Nz`: Grid sizes
- `:tmax`: Requested max time

# Examples
```julia
# Single rank with square grid
results = run_simulation(Nx=20, Ny=20, tmax=0.1)

# Multi-rank (run with: mpiexec -n 4 julia script.jl)
results = run_simulation(Nx=40, Ny=40, tmax=0.2, num_workers=4)

# With visualization and non-square grid
results = run_simulation(Nx=60, Ny=40, tmax=0.1, enable_plots=true)

# Legacy: Use Np for square grids
results = run_simulation(Np=40, tmax=0.1)
```

# See Also
- [`simulation_runner`](@ref): Lower-level simulation function
- [`run_simulation_with_snapshots`](@ref): Simulation with time-series output
"""
function run_simulation(; Np=nothing, Nx=20, Ny=20, Nz=1, tmax=0.1, num_workers=1, verbose=true, save_output=false,
                          enable_plots=false, save_figures=false, output_dir=".",
                          Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5, homogeneous_z=true, debug_output=false)
    # Initialize MPI if not already done
    if !MPI.Initialized()
        MPI.Init()
    end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    # Handle legacy Np parameter (for backward compatibility)
    if Np !== nothing
        Nx = Np
        Ny = Np
    end
    
    # Verify num_workers matches actual MPI size
    if nprocs != num_workers
        if rank == 0
            @warn "num_workers=$num_workers but MPI size=$nprocs. Using $nprocs ranks."
        end
    end
    
    # Setup parameters
    # Domain is [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]
    dx = 1.0 / Nx
    dy = 1.0 / Ny
    dz = 1.0 / Nz
    dtmax = CFL * min(dx, dy, dz)
    nnmax = ceil(Int, tmax / dtmax) + 100000  # Safety margin
    
    # Initial condition parameters (crossing jets)
    rhol = 1.0    # High density in jets
    rhor = 0.01   # Low density background
    T = 1.0
    r110 = 0.0
    r101 = 0.0
    r011 = 0.0
    
    params = (
        Nx = Nx,
        Ny = Ny,
        Nz = Nz,
        tmax = tmax,
        Kn = Kn,
        Ma = Ma,
        flag2D = flag2D,
        CFL = CFL,
        dx = dx,
        dy = dy,
        dz = dz,
        Nmom = 35,
        nnmax = nnmax,
        dtmax = dtmax,
        rhol = rhol,
        rhor = rhor,
        T = T,
        r110 = r110,
        r101 = r101,
        r011 = r011,
        symmetry_check_interval = 10,
        homogeneous_z = homogeneous_z,
        enable_memory_tracking = false,
        debug_output = debug_output
    )
    
    if verbose && rank == 0
        println("="^60)
        println("HyQMOM Simulation")
        println("="^60)
        println("Grid: $(Nx)x$(Ny)x$(Nz)")
        println("MPI ranks: $nprocs")
        println("tmax: $tmax")
        println("Kn: $Kn, Ma: $Ma")
        println("CFL: $CFL, dtmax: $dtmax")
        println("homogeneous_z: $homogeneous_z")
        println("="^60)
    end
    
    # Run simulation
    M_final, final_time, time_steps, grid_out = simulation_runner(params)
    
    if verbose && rank == 0
        println("="^60)
        println("Simulation complete!")
        println("Final time: $final_time")
        println("Time steps: $time_steps")
        println("="^60)
    end
    
    # Generate plots if requested (only on rank 0)
    if enable_plots && rank == 0 && !isnothing(M_final)
        plot_final_results(M_final, grid_out.xm, grid_out.ym, Nx, 35; 
                         save_figures=save_figures, output_dir=output_dir,
                         zm=grid_out.zm, Nz=Nz)
    end
    
    # Package results
    if rank == 0
        return Dict(
            :M => M_final,
            :final_time => final_time,
            :time_steps => time_steps,
            :xm => grid_out.xm,
            :ym => grid_out.ym,
            :zm => grid_out.zm,
            :Nx => Nx,
            :Ny => Ny,
            :Nz => Nz,
            :tmax => tmax
        )
    else
        return Dict(
            :M => nothing,
            :final_time => final_time,
            :time_steps => time_steps,
            :xm => nothing,
            :ym => nothing,
            :zm => nothing,
            :Nx => Nx,
            :Ny => Ny,
            :Nz => Nz,
            :tmax => tmax
        )
    end
end

"""
    run_simulation_with_snapshots(params; snapshot_interval=2)

Run simulation with time-series snapshot collection.

This function extends [`simulation_runner`](@ref) to automatically collect and save
simulation snapshots at regular intervals. Snapshots are streamed to a JLD2 file
for memory efficiency and can be visualized with `interactive_3d_timeseries_streaming`.

# Arguments
- `params`: Simulation parameters (same as [`simulation_runner`](@ref))
- `snapshot_interval`: Save snapshots every N time steps (default: 2)

# Returns
- `filename`: Path to JLD2 snapshot file (rank 0 only)
- `grid`: Grid structure for visualization (rank 0 only)

On non-zero ranks, returns `(nothing, nothing)`.

# Examples
```julia
using HyQMOM, MPI

MPI.Init()

# Set up parameters
params = (
    Nx = 40, Ny = 40, Nz = 40,
    tmax = 0.1, Ma = 1.0, Kn = 1.0,
    CFL = 0.7, flag2D = 0,
    # ... other parameters
)

# Run with snapshot collection
filename, grid = run_simulation_with_snapshots(params; snapshot_interval=5)

# Visualize results (rank 0 only)
if MPI.Comm_rank(MPI.COMM_WORLD) == 0 && filename !== nothing
    interactive_3d_timeseries_streaming(filename, grid, params)
end

MPI.Finalize()
```

# See Also
- [`simulation_runner`](@ref): Core simulation function
- `interactive_3d_timeseries_streaming`: Visualization of snapshot files (requires GLMakie)
"""
function run_simulation_with_snapshots(params; snapshot_interval=2)
    # Add snapshot parameters to the input params
    params_with_snapshots = merge(params, (snapshot_interval=snapshot_interval,))
    
    # Run simulation with snapshot collection
    result = simulation_runner(params_with_snapshots)
    
    # simulation_runner returns different things based on snapshot_interval
    if haskey(params_with_snapshots, :snapshot_interval) && params_with_snapshots.snapshot_interval > 0
        # Snapshot mode: returns (filename, grid) or (nothing, nothing)
        return result
    else
        # Standard mode: returns (M_final, final_time, time_steps, grid)
        # Convert to snapshot-like return format
        M_final, final_time, time_steps, grid = result
        return nothing, grid  # No snapshot file created
    end
end
