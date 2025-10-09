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
  - `Np`: Global grid size (Np × Np)
  - `tmax`: Maximum simulation time
  - `Kn`: Knudsen number
  - `Ma`: Mach number
  - `flag2D`: 2D flag (1 for 2D, 0 for 3D)
  - `CFL`: CFL number
  - `dx`, `dy`: Grid spacing
  - `Nmom`: Number of moments (35)
  - `nnmax`: Maximum number of time steps
  - `dtmax`: Maximum time step size
  - IC parameters: `rhol`, `rhor`, `T`, `r110`, `r101`, `r011`
  - `symmetry_check_interval`: How often to check symmetry
  - `enable_memory_tracking`: Enable memory tracking (not implemented)

# Returns
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
    Np = params.Np
    tmax = params.tmax
    Kn = params.Kn
    Ma = params.Ma
    flag2D = params.flag2D
    CFL = params.CFL
    dx = params.dx
    dy = params.dy
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
    
    # Setup domain decomposition
    decomp = setup_mpi_cartesian_2d(Np, halo, comm)
    nx = decomp.local_size[1]
    ny = decomp.local_size[2]
    
    # Create local grid structure
    xmin, xmax = -0.5, 0.5
    ymin, ymax = -0.5, 0.5
    dx_global = (xmax - xmin) / Np
    dy_global = (ymax - ymin) / Np
    
    # Override params.dx/dy with correct grid spacing
    # (params.dx may have been set incorrectly as 2.0/Np in wrapper)
    dx = dx_global
    dy = dy_global
    
    # Local grid coordinates
    i0i1 = decomp.istart_iend
    j0j1 = decomp.jstart_jend
    
    # Cell edges for local subdomain
    x_local = range(xmin + (i0i1[1]-1)*dx_global, step=dx_global, length=nx+1)
    y_local = range(ymin + (j0j1[1]-1)*dy_global, step=dy_global, length=ny+1)
    
    # Cell centers
    xm_local = collect(x_local[1:end-1]) .+ dx_global/2
    ym_local = collect(y_local[1:end-1]) .+ dy_global/2
    
    grid_local = (dx=dx_global, dy=dy_global,
                  x=collect(x_local), y=collect(y_local),
                  xm=xm_local, ym=ym_local)
    
    # Allocate local arrays with halos
    M = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
    Mnp = similar(M)
    Fx = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
    Fy = zeros(Float64, nx+2*halo, ny+2*halo, Nmom)
    
    # Wave speed arrays (interior only)
    vpxmin = zeros(Float64, nx, ny)
    vpxmax = zeros(Float64, nx, ny)
    vpymin = zeros(Float64, nx, ny)
    vpymax = zeros(Float64, nx, ny)
    v5xmin = zeros(Float64, nx, ny)
    v5xmax = zeros(Float64, nx, ny)
    v5ymin = zeros(Float64, nx, ny)
    v5ymax = zeros(Float64, nx, ny)
    v6xmin = zeros(Float64, nx, ny)
    v6xmax = zeros(Float64, nx, ny)
    v6ymin = zeros(Float64, nx, ny)
    v6ymax = zeros(Float64, nx, ny)
    
    # Build initial conditions locally
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
    Csize = floor(Int, 0.1 * Np)
    Mint = div(Np, 2) + 1
    Maxt = div(Np, 2) + 1 + Csize
    Minb = div(Np, 2) - Csize
    Maxb = div(Np, 2)
    
    # Fill local subdomain with appropriate IC
    for ii in 1:nx
        gi = i0i1[1] + ii - 1  # global i index
        for jj in 1:ny
            gj = j0j1[1] + jj - 1  # global j index
            
            # Default: background
            Mr = Mr_bg
            
            # Bottom jet (moving up-right)
            if gi >= Minb && gi <= Maxb && gj >= Minb && gj <= Maxb
                Mr = Mb
            end
            
            # Top jet (moving down-left) - overwrites if overlapping
            if gi >= Mint && gi <= Maxt && gj >= Mint && gj <= Maxt
                Mr = Mt
            end
            
            M[ii + halo, jj + halo, :] = Mr
        end
    end
    
    # Initial halo exchange
    halo_exchange_2d!(M, decomp, bc)
    
    # Time evolution
    t = 0.0
    nn = 0
    
    if rank == 0
        println("Starting time evolution...")
        println("  Grid: $(Np)×$(Np), Ranks: $(nprocs), Local: $(nx)×$(ny)")
        println("  tmax: $(tmax), CFL: $(CFL), Ma: $(Ma), Kn: $(Kn)")
    end
    
    while t < tmax && nn < nnmax
        nn += 1
        step_start_time = time()
        
        # DEBUG: Track cell (7,13) - M[3]
        debug_cell_i = 7
        debug_cell_j = 13
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            debug_ih = debug_cell_i + halo
            debug_jh = debug_cell_j + halo
            @printf("\n[Step %d] Cell (%d,%d) M[3] tracking:\n", nn, debug_cell_i, debug_cell_j)
            @printf("  Start: M[3] = %.6e\n", M[debug_ih, debug_jh, 3])
        end
        
        # Compute fluxes and wave speeds for interior cells
        for i in 1:nx
            for j in 1:ny
                ih = i + halo
                jh = j + halo
                MOM = M[ih, jh, :]
                
                # Debug label for cell (7,13)
                cell_label_flux = (i == debug_cell_i && j == debug_cell_j) ? " @cell(7,13)" : ""
                
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma, debug_label="[flux_comp1]$(cell_label_flux)")
                v6xmin[i,j], v6xmax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6ymin[i,j], v6ymax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                Mx, My, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma, debug_label="[flux_comp2]$(cell_label_flux)")
                
                Fx[ih, jh, :] = Mx
                Fy[ih, jh, :] = My
                Mnp[ih, jh, :] = Mr
                
                _, v5xmin[i,j], v5xmax[i,j] = closure_and_eigenvalues(Mr[[1,2,3,4,5]])
                _, v5ymin[i,j], v5ymax[i,j] = closure_and_eigenvalues(Mr[[1,6,10,13,15]])
                
                vpxmin[i,j] = min(v5xmin[i,j], v6xmin[i,j])
                vpxmax[i,j] = max(v5xmax[i,j], v6xmax[i,j])
                vpymin[i,j] = min(v5ymin[i,j], v6ymin[i,j])
                vpymax[i,j] = max(v5ymax[i,j], v6ymax[i,j])
            end
        end
        M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]
        
        # DEBUG: After flux computation
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After flux computation: M[3] = %.6e\n", M[debug_ih, debug_jh, 3])
        end
        
        # Exchange M, Fx, Fy from neighbors
        halo_exchange_2d!(M, decomp, bc)
        halo_exchange_2d!(Fx, decomp, bc)
        halo_exchange_2d!(Fy, decomp, bc)
        
        # Compute fluxes and wave speeds in halo cells
        vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext =
            compute_halo_fluxes_and_wavespeeds!(M, Fx, Fy, vpxmin, vpxmax, vpymin, vpymax,
                                               nx, ny, halo, flag2D, Ma)
        
        # Global reduction for time step
        vmax_local = maximum([abs.(vpxmax); abs.(vpxmin); abs.(vpymax); abs.(vpymin)])
        vmax = MPI.Allreduce(vmax_local, max, comm)
        dt = min(CFL*dx/vmax, dtmax)
        dt = min(dt, tmax-t)
        
        # Debug: print vmax if it's unusually large
        if rank == 0 && (vmax > 100.0 || isnan(vmax) || isinf(vmax))
            @printf("  WARNING: vmax = %.6e (unusually large or invalid)\n", vmax)
            # Find which cell has the problem
            for i in 1:nx
                for j in 1:ny
                    if abs(vpxmax[i,j]) > 100.0 || abs(vpxmin[i,j]) > 100.0 ||
                       abs(vpymax[i,j]) > 100.0 || abs(vpymin[i,j]) > 100.0 ||
                       isnan(vpxmax[i,j]) || isnan(vpxmin[i,j]) ||
                       isnan(vpymax[i,j]) || isnan(vpymin[i,j])
                        @printf("    Cell (%d,%d): vpxmin=%.2e, vpxmax=%.2e, vpymin=%.2e, vpymax=%.2e\n",
                               i, j, vpxmin[i,j], vpxmax[i,j], vpymin[i,j], vpymax[i,j])
                        # Print moments at this cell
                        ih = i + halo
                        jh = j + halo
                        @printf("    M[1:5] = [%.6e, %.6e, %.6e, %.6e, %.6e]\n",
                               M[ih,jh,1], M[ih,jh,2], M[ih,jh,3], M[ih,jh,4], M[ih,jh,5])
                    end
                end
            end
        end
        
        t += dt
        
        # X-direction flux update
        Mnpx = apply_flux_update(M, Fx, vpxmin, vpxmax, vpxmin_ext, vpxmax_ext,
                                 nx, ny, halo, dt, dx, decomp, 1)
        
        # DEBUG: After X-flux update
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After X-flux update: Mnpx[3] = %.6e\n", Mnpx[debug_ih, debug_jh, 3])
        end
        
        # Y-direction flux update
        Mnpy = apply_flux_update(M, Fy, vpymin, vpymax, vpymin_ext, vpymax_ext,
                                 nx, ny, halo, dt, dy, decomp, 2)
        
        # DEBUG: After Y-flux update
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After Y-flux update: Mnpy[3] = %.6e\n", Mnpy[debug_ih, debug_jh, 3])
        end
        
        # Combine updates (Strang splitting) - INTERIOR ONLY
        M[halo+1:halo+nx, halo+1:halo+ny, :] =
            Mnpx[halo+1:halo+nx, halo+1:halo+ny, :] +
            Mnpy[halo+1:halo+nx, halo+1:halo+ny, :] -
            M[halo+1:halo+nx, halo+1:halo+ny, :]
        
        # DEBUG: After combining updates
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After combining (Mnpx+Mnpy-M): M[3] = %.6e\n", M[debug_ih, debug_jh, 3])
        end
        
        # Exchange halos before realizability enforcement
        halo_exchange_2d!(M, decomp, bc)
        
        # DEBUG: After halo exchange
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After halo exchange: M[3] = %.6e\n", M[debug_ih, debug_jh, 3])
        end
        
        # Enforce realizability
        for i in 1:nx
            for j in 1:ny
                ih = i + halo
                jh = j + halo
                MOM = M[ih, jh, :]
                
                # Debug label for cell (7,13)
                cell_label = (i == debug_cell_i && j == debug_cell_j) ? " @cell(7,13)" : ""
                
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma, debug_label="[call1]$(cell_label)")
                v6xmin[i,j], v6xmax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6ymin[i,j], v6ymax[i,j], Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                _, _, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma, debug_label="[call2]$(cell_label)")
                
                # DEBUG: Check for large values before assignment
                if rank == 0 && abs(Mr[3]) > 1e6
                    @printf("  ⚠️  Cell (%d,%d): Mr[3] = %.6e before assignment!\n", i, j, Mr[3])
                end
                
                Mnp[ih, jh, :] = Mr
                
                # DEBUG: Check for large values after assignment
                if rank == 0 && abs(Mnp[ih, jh, 3]) > 1e6
                    @printf("  ⚠️  Cell (%d,%d): Mnp[%d,%d,3] = %.6e after assignment!\n", 
                           i, j, ih, jh, Mnp[ih, jh, 3])
                end
            end
        end
        # DEBUG: Check Mnp before bulk assignment
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  Before bulk assignment: Mnp[%d,%d,3] = %.6e\n", debug_ih, debug_jh, Mnp[debug_ih, debug_jh, 3])
            # Check if it's already corrupted
            if abs(Mnp[debug_ih, debug_jh, 3]) > 1e6
                println("  ❌ ALREADY CORRUPTED before bulk assignment!")
                # Scan all cells to find which one corrupted it
                for ii in 1:nx+2*halo
                    for jj in 1:ny+2*halo
                        if abs(Mnp[ii, jj, 3]) > 1e6
                            @printf("    Cell Mnp[%d,%d,3] = %.6e\n", ii, jj, Mnp[ii, jj, 3])
                        end
                    end
                end
            end
        end
        
        M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]
        
        # DEBUG: Check M after bulk assignment
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After bulk assignment: M[%d,%d,3] = %.6e\n", debug_ih, debug_jh, M[debug_ih, debug_jh, 3])
        end
        
        # DEBUG: After realizability enforcement
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After realizability: M[3] = %.6e\n", M[debug_ih, debug_jh, 3])
        end
        
        # Apply BGK collision
        for i in 1:nx
            for j in 1:ny
                ih = i + halo
                jh = j + halo
                MOM = M[ih, jh, :]
                MMC = collision35(MOM, dt, Kn)
                Mnp[ih, jh, :] = MMC
                
            end
        end
        M[halo+1:halo+nx, halo+1:halo+ny, :] = Mnp[halo+1:halo+nx, halo+1:halo+ny, :]
        
        # DEBUG: After collision
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After collision: M[3] = %.6e\n", M[debug_ih, debug_jh, 3])
        end
        
        # Exchange halos for next iteration
        halo_exchange_2d!(M, decomp, bc)
        
        # DEBUG: After final halo exchange
        if rank == 0 && debug_cell_i <= nx && debug_cell_j <= ny
            @printf("  After final halo exchange: M[3] = %.6e\n", M[debug_ih, debug_jh, 3])
        end
        
        # Symmetry check (simplified - just print timing)
        step_time = time() - step_start_time
        local_grid_points = nx * ny
        time_per_point_local = step_time / local_grid_points
        max_time_per_point = MPI.Allreduce(time_per_point_local, max, comm)
        
        # Print timestep info
        if rank == 0
            # Always print for debugging (was: mod(nn, symmetry_check_interval) == 0 || nn == 1)
            @printf("Step %4d: t = %.6f, dt = %.6e, max s/pt = %.6e s\n",
                   nn, t, dt, max_time_per_point)
        end
    end
    
    if rank == 0
        println("Time evolution complete: $(nn) steps, t = $(t)")
    end
    
    # Gather results to rank 0
    M_interior = M[halo+1:halo+nx, halo+1:halo+ny, :]
    
    if rank == 0
        M_final = gather_M(M_interior, i0i1, j0j1, Np, Nmom, comm)
        final_time = t
        time_steps = nn
        
        # Create global grid structure
        grid_out = (x = collect(range(xmin, xmax, length=Np+1)),
                   y = collect(range(ymin, ymax, length=Np+1)),
                   dx = dx_global,
                   dy = dy_global,
                   xm = collect(range(xmin + dx_global/2, step=dx_global, length=Np)),
                   ym = collect(range(ymin + dy_global/2, step=dy_global, length=Np)))
        
        return M_final, final_time, time_steps, grid_out
    else
        # Other ranks send to rank 0
        send_M(M_interior, i0i1, j0j1, 0, comm)
        return nothing, t, nn, nothing
    end
end

"""
    run_simulation(; Np=20, tmax=0.1, num_workers=1, kwargs...)

High-level wrapper for running simulations with sensible defaults.

# Arguments
- `Np`: Global grid size (Np × Np), default 20
- `tmax`: Maximum simulation time, default 0.1
- `num_workers`: Number of MPI ranks (must match mpiexec -n), default 1
- `verbose`: Print progress information, default true
- `save_output`: Save output to file, default false
- `Kn`: Knudsen number, default 1.0
- `Ma`: Mach number, default 0.0
- `flag2D`: 2D flag (1 for 2D, 0 for 3D), default 1
- `CFL`: CFL number, default 0.5

# Returns
Dictionary with:
- `:M`: Final moment field (Np×Np×35) on rank 0, nothing on other ranks
- `:final_time`: Actual final time reached
- `:time_steps`: Number of time steps taken
- `:xm`, `:ym`: Grid coordinates (only on rank 0)
- `:Np`: Grid size
- `:tmax`: Requested max time

# Example
```julia
# Single rank
results = run_simulation(Np=20, tmax=0.1)

# Multi-rank (run with: mpiexec -n 4 julia script.jl)
results = run_simulation(Np=40, tmax=0.2, num_workers=4)
```
"""
function run_simulation(; Np=20, tmax=0.1, num_workers=1, verbose=true, save_output=false,
                          Kn=1.0, Ma=0.0, flag2D=0, CFL=0.5)
    # Initialize MPI if not already done
    if !MPI.Initialized()
        MPI.Init()
    end
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    
    # Verify num_workers matches actual MPI size
    if nprocs != num_workers
        if rank == 0
            @warn "num_workers=$num_workers but MPI size=$nprocs. Using $nprocs ranks."
        end
    end
    
    # Setup parameters
    # Domain is [-0.5, 0.5] × [-0.5, 0.5], so dx = (0.5 - (-0.5))/Np = 1.0/Np
    dx = 1.0 / Np  # FIX: Was incorrectly 2.0/Np
    dy = 1.0 / Np
    dtmax = CFL * min(dx, dy)
    nnmax = ceil(Int, tmax / dtmax) + 100  # Safety margin
    
    # Initial condition parameters (crossing jets)
    rhol = 1.0    # High density in jets
    rhor = 0.01   # Low density background
    T = 1.0
    r110 = 0.0
    r101 = 0.0
    r011 = 0.0
    
    params = (
        Np = Np,
        tmax = tmax,
        Kn = Kn,
        Ma = Ma,
        flag2D = flag2D,
        CFL = CFL,
        dx = dx,
        dy = dy,
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
        enable_memory_tracking = false
    )
    
    if verbose && rank == 0
        println("="^60)
        println("HyQMOM Simulation")
        println("="^60)
        println("Grid: $(Np)×$(Np)")
        println("MPI ranks: $nprocs")
        println("tmax: $tmax")
        println("Kn: $Kn, Ma: $Ma")
        println("CFL: $CFL, dtmax: $dtmax")
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
    
    # Package results
    if rank == 0
        return Dict(
            :M => M_final,
            :final_time => final_time,
            :time_steps => time_steps,
            :xm => grid_out.xm,
            :ym => grid_out.ym,
            :Np => Np,
            :tmax => tmax
        )
    else
        return Dict(
            :M => nothing,
            :final_time => final_time,
            :time_steps => time_steps,
            :xm => nothing,
            :ym => nothing,
            :Np => Np,
            :tmax => tmax
        )
    end
end
