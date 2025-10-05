function [Fx, Fy, vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext] = ...
    compute_halo_fluxes_and_wavespeeds(M, Fx, Fy, vpxmin, vpxmax, vpymin, vpymax, nx, ny, halo, flag2D, Ma)
%COMPUTE_HALO_FLUXES_AND_WAVESPEEDS Compute fluxes and wave speeds in halo cells

% Syntax:
%   [Fx, Fy, vpxmin_ext, vpxmax_ext, vpymin_ext, vpymax_ext] = ...
%       compute_halo_fluxes_and_wavespeeds(M, Fx, Fy, vpxmin, vpxmax, vpymin, vpymax, ...
%                                          nx, ny, halo, flag2D, Ma)

% Description:
%   After halo exchange, the moment data M in halo cells is available from neighbors.
%   This function computes the corresponding fluxes and wave speeds in those halo cells,
%   which are needed for the pas_HLL stencil at processor boundaries.

% Inputs:
%   M       - Moment array with halos (nx+2*halo, ny+2*halo, Nmom)
%   Fx      - X-flux array with halos (nx+2*halo, ny+2*halo, Nmom)
%   Fy      - Y-flux array with halos (nx+2*halo, ny+2*halo, Nmom)
%   vpxmin  - Min x-wave speed, interior only (nx, ny)
%   vpxmax  - Max x-wave speed, interior only (nx, ny)
%   vpymin  - Min y-wave speed, interior only (nx, ny)
%   vpymax  - Max y-wave speed, interior only (nx, ny)
%   nx      - Interior size in x
%   ny      - Interior size in y
%   halo    - Halo width
%   flag2D  - 2D flag for flux closure
%   Ma      - Mach number for flux closure

% Outputs:
%   Fx          - Updated X-flux array with halo values computed
%   Fy          - Updated Y-flux array with halo values computed
%   vpxmin_ext  - Extended min x-wave speed (nx+2*halo, ny)
%   vpxmax_ext  - Extended max x-wave speed (nx+2*halo, ny)
%   vpymin_ext  - Extended min y-wave speed (nx, ny+2*halo)
%   vpymax_ext  - Extended max y-wave speed (nx, ny+2*halo)

    % Create extended wave speed arrays (interior + halos)
    vpxmin_ext = zeros(nx+2*halo, ny);
    vpxmax_ext = zeros(nx+2*halo, ny);
    vpymin_ext = zeros(nx, ny+2*halo);
    vpymax_ext = zeros(nx, ny+2*halo);
    
    % Copy interior wave speeds
    vpxmin_ext(halo+1:halo+nx, :) = vpxmin;
    vpxmax_ext(halo+1:halo+nx, :) = vpxmax;
    vpymin_ext(:, halo+1:halo+ny) = vpymin;
    vpymax_ext(:, halo+1:halo+ny) = vpymax;
    
    % Compute Fx, Fy, and wave speeds in halo cells (they have M data from exchange)
    
    % Left halo (i=1:halo)
    for i = 1:halo
        for j = 1:ny
            jh = j + halo;
            MOM = squeeze(M(i, jh, :));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6x_min, v6x_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [v6y_min, v6y_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            Fx(i, jh, :) = Mx;
            Fy(i, jh, :) = My;
            [~, v5x_min, v5x_max] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
            vpxmin_ext(i, j) = min(v5x_min, v6x_min);
            vpxmax_ext(i, j) = max(v5x_max, v6x_max);
        end
    end
    
    % Right halo (i=halo+nx+1:nx+2*halo)
    for i = halo+nx+1:nx+2*halo
        for j = 1:ny
            jh = j + halo;
            MOM = squeeze(M(i, jh, :));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6x_min, v6x_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [v6y_min, v6y_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            Fx(i, jh, :) = Mx;
            Fy(i, jh, :) = My;
            [~, v5x_min, v5x_max] = closure_and_eigenvalues(Mr([1,2,3,4,5]));
            vpxmin_ext(i, j) = min(v5x_min, v6x_min);
            vpxmax_ext(i, j) = max(v5x_max, v6x_max);
        end
    end
    
    % Bottom halo (j=1:halo)
    for i = 1:nx
        ih = i + halo;
        for j = 1:halo
            MOM = squeeze(M(ih, j, :));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6x_min, v6x_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [v6y_min, v6y_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            Fx(ih, j, :) = Mx;
            Fy(ih, j, :) = My;
            [~, v5y_min, v5y_max] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
            vpymin_ext(i, j) = min(v5y_min, v6y_min);
            vpymax_ext(i, j) = max(v5y_max, v6y_max);
        end
    end
    
    % Top halo (j=halo+ny+1:ny+2*halo)
    for i = 1:nx
        ih = i + halo;
        for j = halo+ny+1:ny+2*halo
            MOM = squeeze(M(ih, j, :));
            [~,~,~,Mr] = Flux_closure35_and_realizable_3D(MOM, flag2D, Ma);
            [v6x_min, v6x_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 1, flag2D, Ma);
            [v6y_min, v6y_max, Mr] = eigenvalues6_hyperbolic_3D(Mr, 2, flag2D, Ma);
            [Mx, My, ~, Mr] = Flux_closure35_and_realizable_3D(Mr, flag2D, Ma);
            Fx(ih, j, :) = Mx;
            Fy(ih, j, :) = My;
            [~, v5y_min, v5y_max] = closure_and_eigenvalues(Mr([1,6,10,13,15]));
            vpymin_ext(i, j) = min(v5y_min, v6y_min);
            vpymax_ext(i, j) = max(v5y_max, v6y_max);
        end
    end
end

