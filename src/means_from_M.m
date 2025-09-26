function [umean, vmean, wmean] = means_from_M(M)
% means_from_M Extract mean velocities from moment vector
%
% Input:
%   M - moment vector (35 elements)
%
% Output:
%   umean, vmean, wmean - mean velocities in x, y, z directions

idx = moment_indices();
M000 = M(idx.m000);
umean = M(idx.m100) / M000;
vmean = M(idx.m010) / M000;
wmean = M(idx.m001) / M000;

end
