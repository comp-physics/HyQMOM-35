function [lam6xa, lam6xb, lam6ya, lam6yb] = compute_jacobian_eigenvalues(M)
% compute_jacobian_eigenvalues Compute 6x6 Jacobian eigenvalues for all grid points
%
% Input:
%   M - moment array (Np x Np x Nmom)
%
% Outputs:
%   lam6xa, lam6xb, lam6ya, lam6yb - eigenvalue arrays (Np x Np x 6)

[Np, ~, ~] = size(M);
lam6xa = zeros(Np,Np,6);
lam6xb = zeros(Np,Np,6);
lam6ya = zeros(Np,Np,6);
lam6yb = zeros(Np,Np,6);

for i = 1:Np
    for j = 1:Np
        M1 = squeeze(M(i,j,:));
        
        % Extract moment components (using direct indexing instead of many scalars)
        % For jacobian6, we need specific moment orderings
        
        % xa: (m000,m010,m020,m030,m040,m100,m110,m120,m130,m200,m210,m220,m300,m310,m400)
        args_xa = [M1(1), M1(6), M1(10), M1(13), M1(15), ...
                   M1(2), M1(7), M1(11), M1(14), ...
                   M1(3), M1(8), M1(12), M1(4), M1(9), M1(5)];
        J6 = jacobian6(args_xa(1), args_xa(2), args_xa(3), args_xa(4), args_xa(5), ...
                       args_xa(6), args_xa(7), args_xa(8), args_xa(9), ...
                       args_xa(10), args_xa(11), args_xa(12), args_xa(13), args_xa(14), args_xa(15));
        lam6xa(i,j,:) = eig(J6);
        
        % xb: (m000,m001,m002,m003,m004,m100,m101,m102,m103,m200,m201,m202,m300,m301,m400)
        args_xb = [M1(1), M1(16), M1(20), M1(23), M1(25), ...
                   M1(2), M1(17), M1(21), M1(24), ...
                   M1(3), M1(18), M1(22), M1(4), M1(19), M1(5)];
        J6 = jacobian6(args_xb(1), args_xb(2), args_xb(3), args_xb(4), args_xb(5), ...
                       args_xb(6), args_xb(7), args_xb(8), args_xb(9), ...
                       args_xb(10), args_xb(11), args_xb(12), args_xb(13), args_xb(14), args_xb(15));
        lam6xb(i,j,:) = eig(J6);
        
        % ya: (m000,m100,m200,m300,m400,m010,m110,m210,m310,m020,m120,m220,m030,m130,m040)
        args_ya = [M1(1), M1(2), M1(3), M1(4), M1(5), ...
                   M1(6), M1(7), M1(8), M1(9), ...
                   M1(10), M1(11), M1(12), M1(13), M1(14), M1(15)];
        J6 = jacobian6(args_ya(1), args_ya(2), args_ya(3), args_ya(4), args_ya(5), ...
                       args_ya(6), args_ya(7), args_ya(8), args_ya(9), ...
                       args_ya(10), args_ya(11), args_ya(12), args_ya(13), args_ya(14), args_ya(15));
        lam6ya(i,j,:) = eig(J6);
        
        % yb: (m000,m001,m002,m003,m004,m010,m011,m012,m013,m020,m021,m022,m030,m031,m040)
        args_yb = [M1(1), M1(16), M1(20), M1(23), M1(25), ...
                   M1(6), M1(26), M1(32), M1(34), ...
                   M1(10), M1(29), M1(35), M1(13), M1(31), M1(15)];
        J6 = jacobian6(args_yb(1), args_yb(2), args_yb(3), args_yb(4), args_yb(5), ...
                       args_yb(6), args_yb(7), args_yb(8), args_yb(9), ...
                       args_yb(10), args_yb(11), args_yb(12), args_yb(13), args_yb(14), args_yb(15));
        lam6yb(i,j,:) = eig(J6);
    end
end

end
