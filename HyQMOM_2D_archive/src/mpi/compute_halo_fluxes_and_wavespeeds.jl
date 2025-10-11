"""
    compute_halo_fluxes_and_wavespeeds!(M, Fx, Fy, vpxmin, vpxmax, vpymin, vpymax,
                                        nx, ny, halo, flag2D, Ma)

Compute fluxes and wave speeds in halo cells.

# Description
After halo exchange, the moment data M in halo cells is available from neighbors.
This function computes the corresponding fluxes and wave speeds in those halo cells,
which are needed for the pas_HLL stencil at processor boundaries.

# Arguments
- `M`: Moment array with halos (nx+2*halo, ny+2*halo, Nmom)
- `Fx`: X-flux array with halos (nx+2*halo, ny+2*halo, Nmom) - modified in-place
- `Fy`: Y-flux array with halos (nx+2*halo, ny+2*halo, Nmom) - modified in-place
- `vpxmin`: Min x-wave speed, interior only (nx, ny)
- `vpxmax`: Max x-wave speed, interior only (nx, ny)
- `vpymin`: Min y-wave speed, interior only (nx, ny)
- `vpymax`: Max y-wave speed, interior only (nx, ny)
- `nx`: Interior size in x
- `ny`: Interior size in y
- `halo`: Halo width
- `flag2D`: 2D flag for flux closure
- `Ma`: Mach number for flux closure

# Returns
- `vpxmin_ext`: Extended min x-wave speed (nx+2*halo, ny)
- `vpxmax_ext`: Extended max x-wave speed (nx+2*halo, ny)
- `vpymin_ext`: Extended min y-wave speed (nx, ny+2*halo)
- `vpymax_ext`: Extended max y-wave speed (nx, ny+2*halo)

# Algorithm
1. Create extended wave speed arrays
2. Copy interior wave speeds
3. Compute fluxes and wave speeds in left halo
4. Compute fluxes and wave speeds in right halo
5. Compute fluxes and wave speeds in bottom halo
6. Compute fluxes and wave speeds in top halo

# Notes
- Modifies Fx and Fy in-place
- Returns extended wave speed arrays
- Uses closure_and_eigenvalues for 5-moment eigenvalues
- Uses eigenvalues6_hyperbolic_3D for 6-moment eigenvalues
"""
function compute_halo_fluxes_and_wavespeeds!(M::Array{T,3}, Fx::Array{T,3}, Fy::Array{T,3},
                                            vpxmin::Matrix{T}, vpxmax::Matrix{T},
                                            vpymin::Matrix{T}, vpymax::Matrix{T},
                                            nx::Int, ny::Int, halo::Int,
                                            flag2D::Int, Ma::Real, decomp=nothing) where T
    # Create extended wave speed arrays (interior + halos)
    vpxmin_ext = zeros(T, nx+2*halo, ny)
    vpxmax_ext = zeros(T, nx+2*halo, ny)
    vpymin_ext = zeros(T, nx, ny+2*halo)
    vpymax_ext = zeros(T, nx, ny+2*halo)
    
    # Copy interior wave speeds
    vpxmin_ext[halo+1:halo+nx, :] = vpxmin
    vpxmax_ext[halo+1:halo+nx, :] = vpxmax
    vpymin_ext[:, halo+1:halo+ny] = vpymin
    vpymax_ext[:, halo+1:halo+ny] = vpymax
    
    # Compute Fx, Fy, and wave speeds in halo cells
    # CRITICAL: Only compute for halos at PHYSICAL boundaries, not MPI boundaries!
    # At MPI boundaries, Fx/Fy have already been received via halo exchange.
    
    has_left_neighbor = (decomp !== nothing && decomp.neighbors.left != -1)
    has_right_neighbor = (decomp !== nothing && decomp.neighbors.right != -1)
    has_down_neighbor = (decomp !== nothing && decomp.neighbors.down != -1)
    has_up_neighbor = (decomp !== nothing && decomp.neighbors.up != -1)
    
    # Left halo (i=1:halo) - only if NO left MPI neighbor
    if !has_left_neighbor
        for i in 1:halo
            for j in 1:ny
                jh = j + halo
                MOM = M[i, jh, :]
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6y_min, v6y_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                Mx, My, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
                Fx[i, jh, :] = Mx
                Fy[i, jh, :] = My
                _, v5x_min, v5x_max = closure_and_eigenvalues(Mr[[1,2,3,4,5]])
                vpxmin_ext[i, j] = min(v5x_min, v6x_min)
                vpxmax_ext[i, j] = max(v5x_max, v6x_max)
            end
        end
    else
        # Copy wave speeds from exchanged flux data (need to recompute from Fx/Fy)
        for i in 1:halo
            for j in 1:ny
                jh = j + halo
                MOM = M[i, jh, :]
                # Don't recompute Fx/Fy - they came from neighbor
                # But still need wave speeds
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                _, v5x_min, v5x_max = closure_and_eigenvalues(Mr[[1,2,3,4,5]])
                vpxmin_ext[i, j] = min(v5x_min, v6x_min)
                vpxmax_ext[i, j] = max(v5x_max, v6x_max)
            end
        end
    end
    
    # Right halo (i=halo+nx+1:nx+2*halo) - only if NO right MPI neighbor
    if !has_right_neighbor
        for i in halo+nx+1:nx+2*halo
            for j in 1:ny
                jh = j + halo
                MOM = M[i, jh, :]
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6y_min, v6y_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                Mx, My, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
                Fx[i, jh, :] = Mx
                Fy[i, jh, :] = My
                _, v5x_min, v5x_max = closure_and_eigenvalues(Mr[[1,2,3,4,5]])
                vpxmin_ext[i, j] = min(v5x_min, v6x_min)
                vpxmax_ext[i, j] = max(v5x_max, v6x_max)
            end
        end
    else
        for i in halo+nx+1:nx+2*halo
            for j in 1:ny
                jh = j + halo
                MOM = M[i, jh, :]
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                _, v5x_min, v5x_max = closure_and_eigenvalues(Mr[[1,2,3,4,5]])
                vpxmin_ext[i, j] = min(v5x_min, v6x_min)
                vpxmax_ext[i, j] = max(v5x_max, v6x_max)
            end
        end
    end
    
    # Bottom halo (j=1:halo) - only if NO down MPI neighbor
    if !has_down_neighbor
        for i in 1:nx
            ih = i + halo
            for j in 1:halo
                MOM = M[ih, j, :]
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6y_min, v6y_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                Mx, My, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
                Fx[ih, j, :] = Mx
                Fy[ih, j, :] = My
                _, v5y_min, v5y_max = closure_and_eigenvalues(Mr[[1,6,10,13,15]])
                vpymin_ext[i, j] = min(v5y_min, v6y_min)
                vpymax_ext[i, j] = max(v5y_max, v6y_max)
            end
        end
    else
        for i in 1:nx
            ih = i + halo
            for j in 1:halo
                MOM = M[ih, j, :]
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6y_min, v6y_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                _, v5y_min, v5y_max = closure_and_eigenvalues(Mr[[1,6,10,13,15]])
                vpymin_ext[i, j] = min(v5y_min, v6y_min)
                vpymax_ext[i, j] = max(v5y_max, v6y_max)
            end
        end
    end
    
    # Top halo (j=halo+ny+1:ny+2*halo) - only if NO up MPI neighbor
    if !has_up_neighbor
        for i in 1:nx
            ih = i + halo
            for j in halo+ny+1:ny+2*halo
                MOM = M[ih, j, :]
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6y_min, v6y_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                Mx, My, _, Mr = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma)
                Fx[ih, j, :] = Mx
                Fy[ih, j, :] = My
                _, v5y_min, v5y_max = closure_and_eigenvalues(Mr[[1,6,10,13,15]])
                vpymin_ext[i, j] = min(v5y_min, v6y_min)
                vpymax_ext[i, j] = max(v5y_max, v6y_max)
            end
        end
    else
        for i in 1:nx
            ih = i + halo
            for j in halo+ny+1:ny+2*halo
                MOM = M[ih, j, :]
                _, _, _, Mr = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma)
                v6x_min, v6x_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma)
                v6y_min, v6y_max, Mr = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma)
                _, v5y_min, v5y_max = closure_and_eigenvalues(Mr[[1,6,10,13,15]])
                vpymin_ext[i, j] = min(v5y_min, v6y_min)
                vpymax_ext[i, j] = max(v5y_max, v6y_max)
            end
        end
    end
    
    return vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext
end
