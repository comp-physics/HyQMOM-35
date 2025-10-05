function [v6min, v6max, lam6a, lam6b] = compute_jacobian_eigenvalues(moments_a, moments_b)
% compute_jacobian_eigenvalues computes 6x6 Jacobian eigenvalues for two moment sets
% This extracts the common pattern from eigenvalues6x and eigenvalues6y
% Inputs:
%   moments_a - 15-element vector of moments for first Jacobian
%   moments_b - 15-element vector of moments for second Jacobian
% Outputs:
%   v6min, v6max - min and max eigenvalues across both Jacobians
%   lam6a, lam6b - full eigenvalue vectors for both Jacobians
% First Jacobian
J6 = jacobian6(moments_a(1), moments_a(2), moments_a(3), moments_a(4), moments_a(5), ...
               moments_a(6), moments_a(7), moments_a(8), moments_a(9), moments_a(10), ...
               moments_a(11), moments_a(12), moments_a(13), moments_a(14), moments_a(15));
lam6a = eig(J6);
lam6ar = sort(real(lam6a));
v6min = lam6ar(1);
v6max = lam6ar(6);

% Second Jacobian
J6 = jacobian6(moments_b(1), moments_b(2), moments_b(3), moments_b(4), moments_b(5), ...
               moments_b(6), moments_b(7), moments_b(8), moments_b(9), moments_b(10), ...
               moments_b(11), moments_b(12), moments_b(13), moments_b(14), moments_b(15));
lam6b = eig(J6);
lam6br = sort(real(lam6b));
v6min = min([v6min lam6br(1)]);
v6max = max([v6max lam6br(6)]);

end

