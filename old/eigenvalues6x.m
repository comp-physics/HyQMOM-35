function [v6min,v6max] = eigenvalues6x(M)
% eigenvalues of flux Jacobian in x direction
%
% M = [M00,M10,M20,M30,M40,M01,M11,M21,M31,M02,M12,M22,M03,M13,M04]
%
m00 = M(1);
m10 = M(2);
m20 = M(3);
m30 = M(4);
m40 = M(5);
m01 = M(6);
m11 = M(7);
m21 = M(8);
m31 = M(9);
m02 = M(10);
m12 = M(11);
m22 = M(12);
m03 = M(13);
m13 = M(14);
m04 = M(15);
J6 = jacobian6(m00,m01,m02,m03,m04,m10,m11,m12,m13,m20,m21,m22,m30,m31,m40);
lam6 = eig(J6);
lam6 = sort(real(lam6));
v6min = lam6(1);
v6max = lam6(6);
end