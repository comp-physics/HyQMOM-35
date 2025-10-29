"""
    check_3d_symmetry(M, nx, ny, nz, halo=1)

Check symmetry for crossing jets configuration.

For symmetric crossing jets centered at domain center, we expect:
- Density: ρ(i,j,k) ≈ ρ(nx+1-i, ny+1-j, nz+1-k)
- Momentum: u(i,j,k) ≈ -u(nx+1-i, ny+1-j, nz+1-k)

# Arguments
- `M`: Moment array with halos (nx+2*halo, ny+2*halo, nz, Nmom)
- `nx`, `ny`, `nz`: Interior grid dimensions
- `halo`: Halo width (default 1)

# Returns
- `max_diff`: Maximum relative difference in symmetric pairs
"""
function check_3d_symmetry(M::Array{Float64,4}, nx::Int, ny::Int, nz::Int, halo::Int=1)
    max_diff_rho = 0.0
    max_diff_u = 0.0
    max_diff_v = 0.0
    max_diff_w = 0.0
    
    # Check interior cells only
    for k in 1:nz
        for i in 1:nx
            for j in 1:ny
                # Symmetric partner indices
                i_sym = nx + 1 - i
                j_sym = ny + 1 - j
                k_sym = nz + 1 - k
                
                # Get moment indices (with halo offset)
                ih = i + halo
                jh = j + halo
                ih_sym = i_sym + halo
                jh_sym = j_sym + halo
                
                # Density (should be symmetric)
                rho = M[ih, jh, k, 1]
                rho_sym = M[ih_sym, jh_sym, k_sym, 1]
                if abs(rho + rho_sym) > 1e-10  # Avoid division by zero
                    diff = abs(rho - rho_sym) / (0.5 * abs(rho + rho_sym))
                    max_diff_rho = max(max_diff_rho, diff)
                end
                
                # Momentum (should be anti-symmetric)
                u = M[ih, jh, k, 2] / rho
                u_sym = M[ih_sym, jh_sym, k_sym, 2] / rho_sym
                if abs(u - u_sym) > 1e-10
                    diff = abs(u + u_sym) / max(abs(u), abs(u_sym), 1e-10)
                    max_diff_u = max(max_diff_u, diff)
                end
                
                v = M[ih, jh, k, 6] / rho
                v_sym = M[ih_sym, jh_sym, k_sym, 6] / rho_sym
                if abs(v - v_sym) > 1e-10
                    diff = abs(v + v_sym) / max(abs(v), abs(v_sym), 1e-10)
                    max_diff_v = max(max_diff_v, diff)
                end
                
                w = M[ih, jh, k, 16] / rho
                w_sym = M[ih_sym, jh_sym, k_sym, 16] / rho_sym
                if abs(w - w_sym) > 1e-10
                    diff = abs(w + w_sym) / max(abs(w), abs(w_sym), 1e-10)
                    max_diff_w = max(max_diff_w, diff)
                end
            end
        end
    end
    
    return max(max_diff_rho, max_diff_u, max_diff_v, max_diff_w)
end

