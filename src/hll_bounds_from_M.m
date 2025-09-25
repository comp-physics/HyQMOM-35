function [xmin, xmax, ymin, ymax] = hll_bounds_from_M(Mr, idx)
% hll_bounds_from_M Compute HLL bounds from moment vector
%
% Inputs:
%   Mr - 35-element moment vector
%   idx - moment indices structure
%
% Outputs:
%   xmin, xmax, ymin, ymax - HLL bounds for x and y directions

% 1-D hyqmom for x-direction eigenvalues
MOM5_x = Mr(idx.x_moments); % m000 m100 m200 m300 m400
[~, xmin, xmax] = closure_and_eigenvalues(MOM5_x);

% 1-D hyqmom for y-direction eigenvalues  
MOM5_y = Mr(idx.y_moments); % m000 m010 m020 m030 m040
[~, ymin, ymax] = closure_and_eigenvalues(MOM5_y);

end
