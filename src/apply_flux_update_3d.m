function Mnp = apply_flux_update_3d(M, F, vpmin, vpmax, vpmin_ext, vpmax_ext, ...
                                     nx, ny, nz, halo, dt, ds, decomp, axis)
%APPLY_FLUX_UPDATE_3D Apply flux update with processor boundary handling (3D physical space)
%
% Syntax:
%   Mnp = apply_flux_update_3d(M, F, vpmin, vpmax, vpmin_ext, vpmax_ext, ...
%                              nx, ny, nz, halo, dt, ds, decomp, axis)
%
% Description:
%   Unified flux update for X, Y, and Z directions using pas_HLL, with 
%   special handling for processor boundaries. At processor boundaries in x-y, 
%   includes one halo cell so pas_HLL can see neighbor data. Z-direction is
%   simpler since there's no MPI decomposition (Pz=1).
%
% Inputs:
%   M          - Moment array with halos (nx+2*halo, ny+2*halo, nz, Nmom)
%   F          - Flux array (nx+2*halo, ny+2*halo, nz, Nmom)
%   vpmin      - Min wave speed, interior (nx, ny, nz) for all directions
%   vpmax      - Max wave speed, interior (nx, ny, nz) for all directions
%   vpmin_ext  - Min wave speed, extended (nx+2*halo, ny, nz) for X, (nx, ny+2*halo, nz) for Y, unused for Z
%   vpmax_ext  - Max wave speed, extended (nx+2*halo, ny, nz) for X, (nx, ny+2*halo, nz) for Y, unused for Z
%   nx         - Interior size in x
%   ny         - Interior size in y
%   nz         - Interior size in z
%   halo       - Halo width (for x-y only, no halos in z)
%   dt         - Time step size
%   ds         - Grid spacing (dx for X, dy for Y, dz for Z)
%   decomp     - Domain decomposition structure with neighbors field
%   axis       - 1 for X-direction, 2 for Y-direction, 3 for Z-direction
%
% Outputs:
%   Mnp - Updated moment array after flux update (nx+2*halo, ny+2*halo, nz, Nmom)

    % Initialize with current state
    Mnp = M;
    
    if axis == 1
        % X-direction flux update (loop over j, k)
        has_left_neighbor = (decomp.neighbors.left ~= -1);
        has_right_neighbor = (decomp.neighbors.right ~= -1);
        
        for k = 1:nz
            for j = 1:ny
                jh = j + halo;
                
                % Determine array extent: include one halo cell at processor boundaries
                if has_left_neighbor && has_right_neighbor
                    i_start = halo;
                    i_end = halo + nx + 1;
                    vp_start = halo;
                    vp_end = halo + nx + 1;
                    apply_bc_left = false;
                    apply_bc_right = false;
                elseif has_left_neighbor
                    i_start = halo;
                    i_end = halo + nx;
                    vp_start = halo;
                    vp_end = halo + nx;
                    apply_bc_left = false;
                    apply_bc_right = true;
                elseif has_right_neighbor
                    i_start = halo + 1;
                    i_end = halo + nx + 1;
                    vp_start = halo + 1;
                    vp_end = halo + nx + 1;
                    apply_bc_left = true;
                    apply_bc_right = false;
                else
                    i_start = halo + 1;
                    i_end = halo + nx;
                    vp_start = 1;
                    vp_end = nx;
                    apply_bc_left = true;
                    apply_bc_right = true;
                end
                
                % Extract array with appropriate extent
                MOM = squeeze(M(i_start:i_end, jh, k, :));
                FX  = squeeze(F(i_start:i_end, jh, k, :));
                
                % Get wave speeds
                if has_left_neighbor || has_right_neighbor
                    vp_min = vpmin_ext(vp_start:vp_end, j, k);
                    vp_max = vpmax_ext(vp_start:vp_end, j, k);
                else
                    vp_min = vpmin(vp_start:vp_end, j, k);
                    vp_max = vpmax(vp_start:vp_end, j, k);
                end
                
                % Call pas_HLL with BC flags
                MNP = pas_HLL(MOM, FX, dt, ds, vp_min, vp_max, apply_bc_left, apply_bc_right);
                
                % Extract interior portion and write back
                if has_left_neighbor && has_right_neighbor
                    Mnp(halo+1:halo+nx, jh, k, :) = MNP(2:end-1, :);
                elseif has_left_neighbor
                    Mnp(halo+1:halo+nx, jh, k, :) = MNP(2:end, :);
                elseif has_right_neighbor
                    Mnp(halo+1:halo+nx, jh, k, :) = MNP(1:end-1, :);
                else
                    Mnp(halo+1:halo+nx, jh, k, :) = MNP;
                end
            end
        end
        
    elseif axis == 2
        % Y-direction flux update (loop over i, k)
        has_down_neighbor = (decomp.neighbors.down ~= -1);
        has_up_neighbor = (decomp.neighbors.up ~= -1);
        
        for k = 1:nz
            for i = 1:nx
                ih = i + halo;
                
                % Determine array extent: include one halo cell at processor boundaries
                if has_down_neighbor && has_up_neighbor
                    j_start = halo;
                    j_end = halo + ny + 1;
                    vp_start = halo;
                    vp_end = halo + ny + 1;
                    apply_bc_down = false;
                    apply_bc_up = false;
                elseif has_down_neighbor
                    j_start = halo;
                    j_end = halo + ny;
                    vp_start = halo;
                    vp_end = halo + ny;
                    apply_bc_down = false;
                    apply_bc_up = true;
                elseif has_up_neighbor
                    j_start = halo + 1;
                    j_end = halo + ny + 1;
                    vp_start = halo + 1;
                    vp_end = halo + ny + 1;
                    apply_bc_down = true;
                    apply_bc_up = false;
                else
                    j_start = halo + 1;
                    j_end = halo + ny;
                    vp_start = 1;
                    vp_end = ny;
                    apply_bc_down = true;
                    apply_bc_up = true;
                end
                
                % Extract array with appropriate extent
                MOM = squeeze(M(ih, j_start:j_end, k, :));
                FY  = squeeze(F(ih, j_start:j_end, k, :));
                
                % Get wave speeds
                if has_down_neighbor || has_up_neighbor
                    vp_min = vpmin_ext(i, vp_start:vp_end, k)';
                    vp_max = vpmax_ext(i, vp_start:vp_end, k)';
                else
                    vp_min = vpmin(i, vp_start:vp_end, k)';
                    vp_max = vpmax(i, vp_start:vp_end, k)';
                end
                
                % Call pas_HLL with BC flags
                MNP = pas_HLL(MOM, FY, dt, ds, vp_min, vp_max, apply_bc_down, apply_bc_up);
                
                % Extract interior portion and write back
                if has_down_neighbor && has_up_neighbor
                    Mnp(ih, halo+1:halo+ny, k, :) = MNP(2:end-1, :);
                elseif has_down_neighbor
                    Mnp(ih, halo+1:halo+ny, k, :) = MNP(2:end, :);
                elseif has_up_neighbor
                    Mnp(ih, halo+1:halo+ny, k, :) = MNP(1:end-1, :);
                else
                    Mnp(ih, halo+1:halo+ny, k, :) = MNP;
                end
            end
        end
        
    else % axis == 3
        % Z-direction flux update (loop over i, j)
        % Simpler: no MPI decomposition in z, always apply BC at both ends
        for j = 1:ny
            jh = j + halo;
            for i = 1:nx
                ih = i + halo;
                
                % Extract 1D array in z-direction
                % Use reshape to ensure 2D even when nz=1 (squeeze can make it 1D)
                MOM = reshape(M(ih, jh, :, :), [nz, size(M, 4)]);  % (nz, Nmom)
                FZ  = reshape(F(ih, jh, :, :), [nz, size(F, 4)]);  % (nz, Nmom)
                
                % Get wave speeds - use reshape to ensure column vector
                vp_min = reshape(vpmin(i, j, :), [nz, 1]);  % (nz, 1)
                vp_max = reshape(vpmax(i, j, :), [nz, 1]);  % (nz, 1)
                
                % Call pas_HLL with BC flags at both ends (no neighbors in z)
                MNP = pas_HLL(MOM, FZ, dt, ds, vp_min, vp_max, true, true);
                
                % Write back - reshape to match array dimensions
                Mnp(ih, jh, :, :) = reshape(MNP, [1, 1, nz, size(MNP, 2)]);
            end
        end
    end
end


