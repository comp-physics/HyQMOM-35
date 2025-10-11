function [conserved, errors] = verify_conservation(M_initial, M_final, tolerance)
%VERIFY_CONSERVATION Check conservation of mass, momentum, and energy
%   [conserved, errors] = verify_conservation(M_initial, M_final, tolerance)
%   
%   Inputs:
%       M_initial - Initial moment array (Nx x Ny x Nmom) or vector (Nmom x 1)
%       M_final   - Final moment array (same size as M_initial)
%       tolerance - Relative tolerance for conservation check
%   
%   Outputs:
%       conserved - Boolean indicating if all quantities are conserved
%       errors    - Struct with relative errors for each conserved quantity

% Handle both grid arrays and single vectors
if ndims(M_initial) == 3
    % Grid array: sum over spatial dimensions
    M_init_sum = squeeze(sum(sum(M_initial, 1), 2));
    M_final_sum = squeeze(sum(sum(M_final, 1), 2));
else
    % Single vector
    M_init_sum = M_initial;
    M_final_sum = M_final;
end

% Extract conserved quantities (indices in 35-moment vector)
% M000 (mass), M100 (x-momentum), M010 (y-momentum), M001 (z-momentum)
conserved_indices = [1, 2, 6, 16];
conserved_names = {'mass', 'momentum_x', 'momentum_y', 'momentum_z'};

errors = struct();
all_conserved = true;

for i = 1:length(conserved_indices)
    idx = conserved_indices(i);
    name = conserved_names{i};
    
    init_val = M_init_sum(idx);
    final_val = M_final_sum(idx);
    
    % Compute relative error
    if abs(init_val) > eps
        rel_error = abs(final_val - init_val) / abs(init_val);
    else
        rel_error = abs(final_val - init_val);
    end
    
    errors.(name) = rel_error;
    
    if rel_error > tolerance
        all_conserved = false;
        fprintf('WARNING: %s not conserved (rel error = %.6e)\n', name, rel_error);
    end
end

conserved = all_conserved;

end
