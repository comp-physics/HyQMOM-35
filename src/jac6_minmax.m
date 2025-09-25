function [vmin, vmax] = jac6_minmax(M, templates)
% jac6_minmax Compute min/max eigenvalues for 6x6 Jacobian using index templates
%
% Inputs:
%   M - 35-element moment vector
%   templates - cell array of index templates for jacobian6 arguments
%
% Outputs:
%   vmin, vmax - minimum and maximum eigenvalues across all templates

vmin = inf;
vmax = -inf;

for t = 1:numel(templates)
    % Extract moments using template indices
    args = num2cell(M(templates{t}));
    
    % Compute Jacobian eigenvalues
    J6 = jacobian6(args{:});
    lam = sort(real(eig(J6)));
    
    % Update min/max
    vmin = min(vmin, lam(1));
    vmax = max(vmax, lam(end));
end

end
