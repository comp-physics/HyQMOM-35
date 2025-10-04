function Mnpx = apply_flux_update_x(M, Fx, vpxmin, vpxmax, vpxmin_ext, vpxmax_ext, ...
                                     nx, ny, halo, dt, dx, decomp)
%APPLY_FLUX_UPDATE_X Apply X-direction flux update with processor boundary handling
%
% Syntax:
%   Mnpx = apply_flux_update_x(M, Fx, vpxmin, vpxmax, vpxmin_ext, vpxmax_ext, ...
%                              nx, ny, halo, dt, dx, decomp)
%
% Description:
%   Applies the X-direction flux update using pas_HLL, with special handling
%   for processor boundaries. At processor boundaries, includes one halo cell
%   so pas_HLL can see neighbor data. At physical boundaries, applies boundary
%   conditions as needed.
%
% Inputs:
%   M           - Moment array with halos (nx+2*halo, ny+2*halo, Nmom)
%   Fx          - X-flux array with halos (nx+2*halo, ny+2*halo, Nmom)
%   vpxmin      - Min x-wave speed, interior (nx, ny)
%   vpxmax      - Max x-wave speed, interior (nx, ny)
%   vpxmin_ext  - Min x-wave speed, extended (nx+2*halo, ny)
%   vpxmax_ext  - Max x-wave speed, extended (nx+2*halo, ny)
%   nx          - Interior size in x
%   ny          - Interior size in y
%   halo        - Halo width
%   dt          - Time step size
%   dx          - Grid spacing in x
%   decomp      - Domain decomposition structure with neighbors field
%
% Outputs:
%   Mnpx - Updated moment array after X-flux update (nx+2*halo, ny+2*halo, Nmom)

    % Initialize with current state
    Mnpx = M;
    
    % Determine boundary types
    has_left_neighbor = (decomp.neighbors.left ~= -1);
    has_right_neighbor = (decomp.neighbors.right ~= -1);
    
    for j = 1:ny
        jh = j + halo;
        
        % Determine array extent: include one halo cell at processor boundaries
        if has_left_neighbor && has_right_neighbor
            % Both processor boundaries
            i_start = halo;
            i_end = halo + nx + 1;
            vp_start = halo;
            vp_end = halo + nx + 1;
            apply_bc_left = false;
            apply_bc_right = false;
        elseif has_left_neighbor
            % Left processor, right physical
            i_start = halo;
            i_end = halo + nx;
            vp_start = halo;
            vp_end = halo + nx;
            apply_bc_left = false;
            apply_bc_right = true;
        elseif has_right_neighbor
            % Left physical, right processor
            i_start = halo + 1;
            i_end = halo + nx + 1;
            vp_start = halo + 1;
            vp_end = halo + nx + 1;
            apply_bc_left = true;
            apply_bc_right = false;
        else
            % Both physical boundaries (1 rank)
            i_start = halo + 1;
            i_end = halo + nx;
            vp_start = 1;
            vp_end = nx;
            apply_bc_left = true;
            apply_bc_right = true;
        end
        
        % Extract array with appropriate extent
        MOM = squeeze(M(i_start:i_end, jh, :));
        FX  = squeeze(Fx(i_start:i_end, jh, :));
        
        % Get wave speeds (use extended array if we included halos, otherwise interior)
        if has_left_neighbor || has_right_neighbor
            vpx_min = vpxmin_ext(vp_start:vp_end, j);
            vpx_max = vpxmax_ext(vp_start:vp_end, j);
        else
            vpx_min = vpxmin(vp_start:vp_end, j);
            vpx_max = vpxmax(vp_start:vp_end, j);
        end
        
        % Call pas_HLL with BC flags
        MNP = pas_HLL(MOM, FX, dt, dx, vpx_min, vpx_max, apply_bc_left, apply_bc_right);
        
        % Extract interior portion of result and write back
        if has_left_neighbor && has_right_neighbor
            % Result is [halo, interior_1:nx, halo], extract middle nx cells
            Mnpx(halo+1:halo+nx, jh, :) = MNP(2:end-1, :);
        elseif has_left_neighbor
            % Result is [halo, interior_1:nx], extract last nx cells
            Mnpx(halo+1:halo+nx, jh, :) = MNP(2:end, :);
        elseif has_right_neighbor
            % Result is [interior_1:nx, halo], extract first nx cells
            Mnpx(halo+1:halo+nx, jh, :) = MNP(1:end-1, :);
        else
            % Result is [interior_1:nx], use as is
            Mnpx(halo+1:halo+nx, jh, :) = MNP;
        end
    end
end

