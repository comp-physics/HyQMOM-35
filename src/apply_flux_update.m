function Mnp = apply_flux_update(M, F, vpmin, vpmax, vpmin_ext, vpmax_ext, ...
                                  nx, ny, halo, dt, ds, decomp, axis)
%APPLY_FLUX_UPDATE Apply flux update with processor boundary handling

% Syntax:
%   Mnp = apply_flux_update(M, F, vpmin, vpmax, vpmin_ext, vpmax_ext, ...
%                           nx, ny, halo, dt, ds, decomp, axis)

% Description:
%   Unified flux update for both X and Y directions using pas_HLL, with 
%   special handling for processor boundaries. At processor boundaries, 
%   includes one halo cell so pas_HLL can see neighbor data.

% Inputs:
%   M          - Moment array with halos (nx+2*halo, ny+2*halo, Nmom)
%   F          - Flux array with halos (nx+2*halo, ny+2*halo, Nmom)
%   vpmin      - Min wave speed, interior (nx, ny) for X or (nx, ny) for Y
%   vpmax      - Max wave speed, interior (nx, ny) for X or (nx, ny) for Y
%   vpmin_ext  - Min wave speed, extended (nx+2*halo, ny) for X or (nx, ny+2*halo) for Y
%   vpmax_ext  - Max wave speed, extended (nx+2*halo, ny) for X or (nx, ny+2*halo) for Y
%   nx         - Interior size in x
%   ny         - Interior size in y
%   halo       - Halo width
%   dt         - Time step size
%   ds         - Grid spacing (dx for X, dy for Y)
%   decomp     - Domain decomposition structure with neighbors field
%   axis       - 1 for X-direction, 2 for Y-direction

% Outputs:
%   Mnp - Updated moment array after flux update (nx+2*halo, ny+2*halo, Nmom)

    % Initialize with current state
    Mnp = M;
    
    if axis == 1
        % X-direction flux update
        has_left_neighbor = (decomp.neighbors.left ~= -1);
        has_right_neighbor = (decomp.neighbors.right ~= -1);
        
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
            MOM = squeeze(M(i_start:i_end, jh, :));
            FX  = squeeze(F(i_start:i_end, jh, :));
            
            % Get wave speeds
            if has_left_neighbor || has_right_neighbor
                vp_min = vpmin_ext(vp_start:vp_end, j);
                vp_max = vpmax_ext(vp_start:vp_end, j);
            else
                vp_min = vpmin(vp_start:vp_end, j);
                vp_max = vpmax(vp_start:vp_end, j);
            end
            
            % Call pas_HLL with BC flags
            MNP = pas_HLL(MOM, FX, dt, ds, vp_min, vp_max, apply_bc_left, apply_bc_right);
            
            % Extract interior portion and write back
            if has_left_neighbor && has_right_neighbor
                Mnp(halo+1:halo+nx, jh, :) = MNP(2:end-1, :);
            elseif has_left_neighbor
                Mnp(halo+1:halo+nx, jh, :) = MNP(2:end, :);
            elseif has_right_neighbor
                Mnp(halo+1:halo+nx, jh, :) = MNP(1:end-1, :);
            else
                Mnp(halo+1:halo+nx, jh, :) = MNP;
            end
        end
        
    else % axis == 2
        % Y-direction flux update
        has_down_neighbor = (decomp.neighbors.down ~= -1);
        has_up_neighbor = (decomp.neighbors.up ~= -1);
        
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
            MOM = squeeze(M(ih, j_start:j_end, :));
            FY  = squeeze(F(ih, j_start:j_end, :));
            
            % Get wave speeds
            if has_down_neighbor || has_up_neighbor
                vp_min = vpmin_ext(i, vp_start:vp_end)';
                vp_max = vpmax_ext(i, vp_start:vp_end)';
            else
                vp_min = vpmin(i, vp_start:vp_end)';
                vp_max = vpmax(i, vp_start:vp_end)';
            end
            
            % Call pas_HLL with BC flags
            MNP = pas_HLL(MOM, FY, dt, ds, vp_min, vp_max, apply_bc_down, apply_bc_up);
            
            % Extract interior portion and write back
            if has_down_neighbor && has_up_neighbor
                Mnp(ih, halo+1:halo+ny, :) = MNP(2:end-1, :);
            elseif has_down_neighbor
                Mnp(ih, halo+1:halo+ny, :) = MNP(2:end, :);
            elseif has_up_neighbor
                Mnp(ih, halo+1:halo+ny, :) = MNP(1:end-1, :);
            else
                Mnp(ih, halo+1:halo+ny, :) = MNP;
            end
        end
    end
end

