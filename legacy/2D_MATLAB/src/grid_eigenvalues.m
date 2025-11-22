function eig_data = grid_eigenvalues(M, Np, Nmom)
%GRID_EIGENVALUES Compute all eigenvalues for entire grid in one pass
%   eig_data = grid_eigenvalues(M, Np, Nmom)
%   Inputs:
%       M    - Np x Np x Nmom moment array
%       Np   - Grid size
%       Nmom - Number of moments (35)
%   Outputs:
%       eig_data - Struct containing:
%           .lam6x  - Np x Np x 6 eigenvalues for UV plane
%           .lam6y  - Np x Np x 6 eigenvalues for VU plane (Y-direction)
%           .lam6z  - Np x Np x 6 eigenvalues for UW plane
%           .lam6w  - Np x Np x 6 eigenvalues for VW plane
%           .v6xmin, .v6xmax - Min/max envelope for X-direction
%           .v6ymin, .v6ymax - Min/max envelope for Y-direction
%   This function eliminates duplicate eigenvalue computations in:
%   - simulation_plots.m (hyperbolicity plots, eigenvalue plots)
%   - main.m (if eigenvalues are needed for analysis)
%   Example:
%       eig = grid_eigenvalues(M_final, Np, 35);
%       % Use eig.lam6x, eig.v6xmin, etc. directly in plots
%   See also: axis_moment_slice, compute_jacobian_eigenvalues
    % Pre-allocate arrays
    eig_data = struct();
    eig_data.lam6x = zeros(Np, Np, 6);  % UV plane (X-direction primary)
    eig_data.lam6y = zeros(Np, Np, 6);  % VU plane (Y-direction primary)
    eig_data.lam6z = zeros(Np, Np, 6);  % UW plane (X-direction secondary)
    eig_data.lam6w = zeros(Np, Np, 6);  % VW plane (Y-direction secondary)
    
    % Compute eigenvalues at each grid point
    for i = 1:Np
        for j = 1:Np
            % Extract moment vector at this grid point
            M1 = zeros(Nmom, 1);
            for k = 1:Nmom
                M1(k) = M(i, j, k);
            end
            
            % Compute eigenvalues for all four planes using axis_moment_slice
            % UV plane (X-direction, primary)
            moments_uv = axis_moment_slice(M1, 1);
            J6 = jacobian6(moments_uv(1), moments_uv(2), moments_uv(3), moments_uv(4), moments_uv(5), ...
                          moments_uv(6), moments_uv(7), moments_uv(8), moments_uv(9), moments_uv(10), ...
                          moments_uv(11), moments_uv(12), moments_uv(13), moments_uv(14), moments_uv(15));
            eig_data.lam6x(i, j, :) = sort(real(eig(J6)));
            
            % VU plane (Y-direction, primary)
            moments_vu = axis_moment_slice(M1, 2);
            J6 = jacobian6(moments_vu(1), moments_vu(2), moments_vu(3), moments_vu(4), moments_vu(5), ...
                          moments_vu(6), moments_vu(7), moments_vu(8), moments_vu(9), moments_vu(10), ...
                          moments_vu(11), moments_vu(12), moments_vu(13), moments_vu(14), moments_vu(15));
            eig_data.lam6y(i, j, :) = sort(real(eig(J6)));
            
            % UW plane (X-direction, secondary)
            moments_uw = axis_moment_slice(M1, 3);
            J6 = jacobian6(moments_uw(1), moments_uw(2), moments_uw(3), moments_uw(4), moments_uw(5), ...
                          moments_uw(6), moments_uw(7), moments_uw(8), moments_uw(9), moments_uw(10), ...
                          moments_uw(11), moments_uw(12), moments_uw(13), moments_uw(14), moments_uw(15));
            eig_data.lam6z(i, j, :) = sort(real(eig(J6)));
            
            % VW plane (Y-direction, secondary)
            moments_vw = axis_moment_slice(M1, 4);
            J6 = jacobian6(moments_vw(1), moments_vw(2), moments_vw(3), moments_vw(4), moments_vw(5), ...
                          moments_vw(6), moments_vw(7), moments_vw(8), moments_vw(9), moments_vw(10), ...
                          moments_vw(11), moments_vw(12), moments_vw(13), moments_vw(14), moments_vw(15));
            eig_data.lam6w(i, j, :) = sort(real(eig(J6)));
        end
    end
    
    % Compute min/max envelopes for each direction
    eig_data.v6xmin = min(eig_data.lam6x, [], 3);
    eig_data.v6xmax = max(eig_data.lam6x, [], 3);
    eig_data.v6ymin = min(eig_data.lam6y, [], 3);
    eig_data.v6ymax = max(eig_data.lam6y, [], 3);
end

