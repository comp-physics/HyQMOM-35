% 3rd leading principal minor
clear
clc
close all

syms S110 S101 S011 S300 S210 S201 S120 S111 S102 S030 S021 S012 S003 S310 S301 S211 S130 S121 S112 S103 S031 S013 real
syms S400 S220 S202 S040 S022 S004 positive

% S101=0;
% S011=0;
% S201=0;
% S111=0;
% S102=0;
% S021=0;
% S012=0;

dDel1 = -(S011^2 - 2*S011*S101*S110 + S101^2 + S110^2 - 1);

[E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);

% % find lower bounds for S220,S202,S002
% E1_2 = simplify(E1(2,2));
% S220a = solve(E1_2,S220);
% E1_3 = simplify(E1(3,3));
% S202a = solve(E1_3,S202);
% E1_5 = simplify(E1(5,5));
% S022a = solve(E1_5,S022);
% 
% E2_2 = simplify(E2(2,2));
% S202b = solve(E2_2,S202);
% E2_3 = simplify(E2(3,3));
% S220b = solve(E2_3,S220);
% E2_5 = simplify(E2(5,5));
% S022b = solve(E2_5,S022);
% 
% E3_2 = simplify(E3(2,2));
% S220c = solve(E3_2,S220);
% E3_3 = simplify(E3(3,3));
% S022c = solve(E3_3,S022);
% E3_5 = simplify(E3(5,5));
% S202c = solve(E3_5,S202);
% 
% E4_2 = simplify(E4(2,2));
% S202d = solve(E4_2,S202);
% E4_3 = simplify(E4(3,3));
% S022d = solve(E4_3,S022);
% E4_5 = simplify(E4(5,5));
% S220d = solve(E4_5,S220);
% 
% E5_2 = simplify(E5(2,2));
% S022e = solve(E5_2,S022)
% E5_3 = simplify(E5(3,3));
% S220e = solve(E5_3,S220)
% E5_5 = simplify(E5(5,5));
% S202e = solve(E5_5,S202)
% 
% E6_2 = simplify(E6(2,2));
% S022f = solve(E6_2,S022);
% E6_3 = simplify(E6(3,3));
% S202f = solve(E6_3,S202);
% E6_5 = simplify(E6(5,5));
% S220f = solve(E6_5,S220);
% 
% simplify(S202a-S202f);
% 
% S22 = [ S220a, S202a, S022a];
% %matlabFunction(S22,'File','lower_bound_S220')

% S101=0;
% S011=0;
% S201=0;
% S111=0;
% S102=0;
% S021=0;
% S012=0;
% S301=0;
% S211=0;
% S121=0;
% S103=0;
% S031=0;
% S013=0;
% S202=1;
% S022=1;
% 
% 
% S22 = lower_bound_S220(S011,S012,S021,S101,S102,S110,S111,S120,S201,S210)


S310a = solve(E1(1,2),S310);
S301a = solve(E1(1,3),S301);
S220a = solve(E1(1,4),S220);
S211a = solve(E1(1,5),S211)
S202a = solve(E1(1,6),S202);

S211bis = solve(E1(2,3),S211)

% S301b = solve(E2(1,2),S301);
% S310b = solve(E2(1,3),S310);
% S202b = solve(E2(1,4),S202);
% S211b = solve(E2(1,5),S211);
% S220b = solve(E2(1,6),S220);
% 
% S130c = solve(E3(1,2),S130);
% S031c = solve(E3(1,3),S031);
% S220c = solve(E3(1,4),S220);
% S121c = solve(E3(1,5),S121);
% S022c = solve(E3(1,6),S022);
% 
% S103d = solve(E4(1,2),S103);
% S013d = solve(E4(1,3),S013);
% S202d = solve(E4(1,4),S202);
% S112d = solve(E4(1,5),S112);
% S022d = solve(E4(1,6),S022);
% 
% S031e = solve(E5(1,2),S031);
% S130e = solve(E5(1,3),S130);
% S022e = solve(E5(1,4),S022);
% S121e = solve(E5(1,5),S121);
% S220e = solve(E5(1,6),S220);
% 
% S013f = solve(E6(1,2),S013);
% S103f = solve(E6(1,3),S103);
% S022f = solve(E6(1,4),S022);
% S112f = solve(E6(1,5),S112);
% S202f = solve(E6(1,6),S202);


simplify(S211a-S211bis)



% M1 = [ S310a, S301a, S130c, S103d, S031c, S013d, S220a, S202a, S022c, S211a, S121c, S112d];
% matlabFunction(M1,'File','bound_minor1')

% S101=0;
% S011=0;
% S201=0;
% S111=0;
% S102=0;
% S021=0;
% S012=0;
% % S301=0;
% % S211=0;
% % S121=0;
% % S103=0;
% % S031=0;
% % S013=0;
% % S202=1;
% % S022=1;
% % 
% % 
% % M1 = bound_minor1(S003,S011,S012,S021,S030,S101,S102,S110,S111,S120,S201,S210,S300);
% % M1 = simplify(M1)
% 
% %simplify(E1(1,4)-E1(2,2))
% 
% [E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
%                       S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
% 
% simplify(E1(1,4)-E1(2,2))

% S220a = solve(E1(1,4),S220)
% S202a = solve(E1(1,6),S202)
% S022a = solve(E4(1,6),S022)
% 
% 
% S22 = [ S220a, S202a, S022a];
% matlabFunction(S22,'File','upper_bound_S220')