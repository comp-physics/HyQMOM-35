function Mnpy = apply_flux_update_y(M, Fy, vpymin, vpymax, vpymin_ext, vpymax_ext, ...
                                     nx, ny, halo, dt, dy, decomp)
%APPLY_FLUX_UPDATE_Y Apply Y-direction flux update with processor boundary handling
%
% Syntax:
%   Mnpy = apply_flux_update_y(M, Fy, vpymin, vpymax, vpymin_ext, vpymax_ext, ...
%                              nx, ny, halo, dt, dy, decomp)
%
% Description:
%   Applies the Y-direction flux update using pas_HLL, with special handling
%   for processor boundaries. At processor boundaries, includes one halo cell
%   so pas_HLL can see neighbor data. At physical boundaries, applies boundary
%   conditions as needed.
%
% Inputs:
%   M           - Moment array with halos (nx+2*halo, ny+2*halo, Nmom)
%   Fy          - Y-flux array with halos (nx+2*halo, ny+2*halo, Nmom)
%   vpymin      - Min y-wave speed, interior (nx, ny)
%   vpymax      - Max y-wave speed, interior (nx, ny)
%   vpymin_ext  - Min y-wave speed, extended (nx, ny+2*halo)
%   vpymax_ext  - Max y-wave speed, extended (nx, ny+2*halo)
%   nx          - Interior size in x
%   ny          - Interior size in y
%   halo        - Halo width
%   dt          - Time step size
%   dy          - Grid spacing in y
%   decomp      - Domain decomposition structure with neighbors field
%
% Outputs:
%   Mnpy - Updated moment array after Y-flux update (nx+2*halo, ny+2*halo, Nmom)

    % Initialize with current state
    Mnpy = M;
    
    % Determine boundary types
    has_down_neighbor = (decomp.neighbors.down ~= -1);
    has_up_neighbor = (decomp.neighbors.up ~= -1);
    
    for i = 1:nx
        ih = i + halo;
        
        % Determine array extent: include one halo cell at processor boundaries
        if has_down_neighbor && has_up_neighbor
            % Both processor boundaries: include both halos
            j_start = halo;
            j_end = halo + ny + 1;
            vp_start = halo;
            vp_end = halo + ny + 1;
            apply_bc_down = false;
            apply_bc_up = false;
        elseif has_down_neighbor
            % Bottom processor, top physical: include bottom halo only
            j_start = halo;
            j_end = halo + ny;
            vp_start = halo;
            vp_end = halo + ny;
            apply_bc_down = false;
            apply_bc_up = true;
        elseif has_up_neighbor
            % Bottom physical, top processor: include top halo only
            j_start = halo + 1;
            j_end = halo + ny + 1;
            vp_start = halo + 1;
            vp_end = halo + ny + 1;
            apply_bc_down = true;
            apply_bc_up = false;
        else
            % Both physical boundaries (1 rank): interior only
            j_start = halo + 1;
            j_end = halo + ny;
            vp_start = 1;
            vp_end = ny;
            apply_bc_down = true;
            apply_bc_up = true;
        end
        
        % Extract array with appropriate extent
        MOM = squeeze(M(ih, j_start:j_end, :));
        FY  = squeeze(Fy(ih, j_start:j_end, :));
        
        % Get wave speeds (use extended array if we included halos, otherwise interior)
        if has_down_neighbor || has_up_neighbor
            vpy_min = vpymin_ext(i, vp_start:vp_end)';
            vpy_max = vpymax_ext(i, vp_start:vp_end)';
        else
            vpy_min = vpymin(i, vp_start:vp_end)';
            vpy_max = vpymax(i, vp_start:vp_end)';
        end
        
        % Call pas_HLL with BC flags
        MNP = pas_HLL(MOM, FY, dt, dy, vpy_min, vpy_max, apply_bc_down, apply_bc_up);
        
        % Extract interior portion of result and write back
        if has_down_neighbor && has_up_neighbor
            % Result is [halo, interior_1:ny, halo], extract middle ny cells
            Mnpy(ih, halo+1:halo+ny, :) = MNP(2:end-1, :);
        elseif has_down_neighbor
            % Result is [halo, interior_1:ny], extract last ny cells
            Mnpy(ih, halo+1:halo+ny, :) = MNP(2:end, :);
        elseif has_up_neighbor
            % Result is [interior_1:ny, halo], extract first ny cells
            Mnpy(ih, halo+1:halo+ny, :) = MNP(1:end-1, :);
        else
            % Result is [interior_1:ny], use as is
            Mnpy(ih, halo+1:halo+ny, :) = MNP;
        end
    end
end

