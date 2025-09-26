function [v6min, v6max, lam6a, lam6b] = compute_eigenvalues_from_moments(M, direction)
% compute_eigenvalues_from_moments Compute min/max eigenvalues using templates
%
% Inputs:
%   M - moment vector (35 elements)
%   direction - 'x' or 'y'
%
% Outputs:
%   v6min, v6max - minimum and maximum eigenvalues
%   lam6a, lam6b - eigenvalue arrays for both Jacobian matrices

templates = jacobian_templates(direction);
[v6min, v6max] = jac6_minmax(M, templates);

% Get individual eigenvalue arrays if requested
if nargout > 2
    args_a = num2cell(M(templates{1}));
    J6a = jacobian6(args_a{:});
    lam6a = eig(J6a);
    
    args_b = num2cell(M(templates{2}));
    J6b = jacobian6(args_b{:});
    lam6b = eig(J6b);
end

end
