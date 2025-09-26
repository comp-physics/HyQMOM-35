function [sC200, sC020, sC002] = sC_from_C4(C4)
% sC_from_C4 Extract square roots of second-order central moments
%
% Input:
%   C4 - central moment array (35 elements)
%
% Output:
%   sC200, sC020, sC002 - square roots of diagonal second-order central moments

idx = moment_indices();
sC200 = sqrt(max(eps, C4(idx.C200)));
sC020 = sqrt(max(eps, C4(idx.C020)));
sC002 = sqrt(max(eps, C4(idx.C002)));

end
