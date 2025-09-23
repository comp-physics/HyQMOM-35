% independint formulas for determinant
clear
clc

% syms S110 S101 S011 S300 S210 S201 S120 S111 S102 S030 S021 S012 S003 S310 S301 S211 S130 S121 S112 S103 S031 S013 real
% syms S400 S220 S202 S040 S022 S004 positive
% 
% dDel1 = -(S011^2 - 2*S011*S101*S110 + S101^2 + S110^2 - 1);
% 
% [E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
%                       S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
% %
% 
% E1 = simplify(E1);
% 
% dE1 = det(E1);
% 
% diff(dE1,S202)
% diff(dE1,S022)
% diff(diff(dE1,S202),S202)
% C1 = diff(diff(diff(dE1,S220),S220),S220);
% C1 = collect(simplify(C1/6),S202)

% ex = simplify(E1(1,4)-E1(2,2));
% ey = simplify(E1(1,6)-E1(3,3));
% ez = simplify(E1(4,6)-E1(5,5));

syms X Y Z e11 e44 e66 positive
syms e12 e13 e15 e23 e24 e25 e26 e34 e35 e36 e45 e56 ex ey ez real

D = [e11  e12 e13 X+ex e15 Y+ey;...
     e12  X   e23 e24  e25 e26 ;...
     e13  e23 Y   e34  e35 e36 ;...
     X+ex e24 e34 e44  e45 Z+ez;...
     e15  e25 e35 e45  Z   e56 ;...
     Y+ey e26 e36 Z+ez e56 e66]
%dD4 = eig(D(1:4,1:4))
%dD5 = eig(D(1:5,1:5))
dD6 = hermiteForm(D,X)
%dD3 = collect(dD3,[Z^3 Z^2 Z])

%R = solve(dD3,Z,MaxDegree=3);

%matlabFunction(R,'File','rootsZ_3D')
%matlabFunction(dD3,'File','detD3_3D')
