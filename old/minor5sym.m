% 3 independint formulas for 4th leading principal minors (all permutations)
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
% E2 = simplify(E2);
% E3 = simplify(E3);
% E4 = simplify(E4);
% E5 = simplify(E5);
% E6 = simplify(E6);
% 
% d4E1 = det(E1(1:4,1:4));
% d4E3 = det(E3(1:4,1:4));

% diff(d4E1,S202)
% diff(d4E1,S022)
% diff(diff(d4E1,S202),S202)
% C1 = diff(diff(diff(d4E1,S220),S220),S220);
% C1 = collect(simplify(C1/6),S202)

% diff(d4E3,S022)
% diff(d4E3,S202)
% diff(diff(d4E3,S022),S022)
% C3 = diff(diff(diff(d4E3,S220),S220),S220);

syms X Y Z e11 e44 positive
syms e12 e13 e15 e23 e24 e25 e34 e35 e45 ex real

% ex = simplify(E1(1,4)-E1(2,2));
% Y = E1(2,2);

D1 = [e11  e12 e13 X+ex e15; ...
      e12  X   e23 e24  e25; ...
      e13  e23 Y   e34  e35; ...
      X+ex e24 e34 e44  e45; ...
      e15  e25 e35 e45  Z];
dD1 = det(D1);
dD1 = collect(dD1,X)

R = solve(dD1,X,MaxDegree=3);
%R = double(R);
matlabFunction(R,'File','rootsR_X_YZ')
