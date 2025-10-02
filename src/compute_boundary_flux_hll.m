function flux = compute_boundary_flux_hll(M_left, M_right, F_left, F_right, vpmin_left, vpmax_left, vpmin_right, vpmax_right)
% compute_boundary_flux_hll - Compute HLL flux at a processor boundary
%
% Inputs:
%   M_left, M_right - Moment vectors on left and right sides of interface (row vectors or columns)
%   F_left, F_right - Flux vectors on left and right sides
%   vpmin_left, vpmax_left - Wave speeds on left side (scalars)
%   vpmin_right, vpmax_right - Wave speeds on right side (scalars)
%
% Output:
%   flux - HLL flux at the interface (same shape as M_left)

    % Ensure row vectors for computation
    M_left = M_left(:)';
    M_right = M_right(:)';
    F_left = F_left(:)';
    F_right = F_right(:)';
    
    % Compute wave speeds at interface
    lleft = min(vpmin_left, vpmin_right);
    lright = max(vpmax_left, vpmax_right);
    
    % Compute Wstar at interface
    if abs(lleft - lright) > 1.d-10
        Wstar = (lleft*M_left - lright*M_right)/(lleft - lright) - (F_left - F_right)/(lleft - lright);
    else
        Wstar = zeros(size(M_left));
    end
    
    % Compute HLL flux (same as flux_HLL for a single interface)
    flux = 0.5*(F_left + F_right) - 0.5*( abs(lleft)*(Wstar - M_left) - abs(lright)*(Wstar - M_right) );
    
    % Return as column vector
    flux = flux(:);
    
end

