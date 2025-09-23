% realizability regions for S220, S202, S002
clear
clc
close all

% select univariate moments
S300=0;
S030=0;
S003=0;
S400=3;
S040=3;
S004=3;
%
H200=S400-S300^2-1;
H020=S040-S030^2-1;
H002=S004-S003^2-1;
% initial guesses based on absolute max
s220max = 1 + sqrt((S400-1)*(S040-1));
s202max = 1 + sqrt((S400-1)*(S004-1));
s022max = 1 + sqrt((S040-1)*(S004-1));
%
S220 = 0.9*s220max
S202 = 0.2*s202max
S022 = 0.2*s022max
%
% select covariances
S110= 0.3;
S101= 0.;
S011= 0.;
[S110,S101,S011] = realizability_S2(S110,S101,S011)
%
Del1 = [1 0    0    0;...
        0 1    S110 S101;...
        0 S110 1    S011;...
        0 S101 S011 1];
%
% select third-order moments
S210=S110*S300;
S201=S101*S300;
S120=S110*S030;
S102=S101*S003;
S021=S011*S030;
S012=S011*S003;
S111 = 0;
%
% 3-order moments:
% etc.
[S210,S201] = realizability_S210(S110,S101,S011,S300,S210,S201,H200);
% check and correct realizability of S120, S021
[S120,S021] = realizability_S210(S110,S011,S101,S030,S120,S021,H020);
% check and correct realizability of S102, S012
[S102,S012] = realizability_S210(S101,S011,S110,S003,S102,S012,H002);
% check and correct realizability of S111
[S111] = realizability_S111(S110,S101,S011,S210,S201,S120,S021,S102,S012,S111);
% find D matrix
B = [1    S300 S210 S201;...
     S110 S210 S120 S111;...
     S101 S201 S111 S102;...
     1    S120 S030 S021;...
     S011 S111 S021 S012;...
     1    S102 S012 S003];
D = B/Del1*B';
%
d11 = D(1,1);
d22 = D(2,2);
d33 = D(3,3);
d44 = D(4,4);
d55 = D(5,5);
d66 = D(6,6);
d12 = D(1,2);
d13 = D(1,3);
d14 = D(1,4);
d15 = D(1,5);
d16 = D(1,6);
d23 = D(2,3);
d24 = D(2,4);
d25 = D(2,5);
d26 = D(2,6);
d34 = D(3,4);
d35 = D(3,5);
d36 = D(3,6);
d45 = D(4,5);
d46 = D(4,6);
d56 = D(5,6);

% select fourth-order moments
S310=S110*S400;
S301=S101*S400;
S130=S110*S040;
S103=S101*S004;
S031=S011*S040;
S013=S011*S004;
%
[S310,S220] = realizability_S310(S110,S101,S011,S300,S210,S201,S120,S111,S310,S220,H200);
[S130,~]    = realizability_S310(S110,S011,S101,S030,S120,S021,S210,S111,S130,S220,H020);
[S301,S202] = realizability_S310(S101,S110,S011,S300,S201,S210,S102,S111,S301,S202,H200);
[S103,~]    = realizability_S310(S101,S011,S110,S003,S102,S012,S201,S111,S103,S202,H002);
[S031,S022] = realizability_S310(S011,S110,S101,S030,S021,S120,S012,S111,S031,S022,H020);
[S013,~]    = realizability_S310(S011,S101,S110,S003,S012,S102,S021,S111,S013,S022,H002);
%
S211=0;
S121=0;
S112=0;
[E1,~,E3,~,~,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
e11a = E1(1,1);
e12a = E1(1,2);
e13a = E1(1,3);
d22a = S220-E1(2,2);
d23a = S211-E1(2,3);
d33a = S202-E1(3,3);
S220a = d22a+e13a^2/e11a;
S202a = d33a+e12a^2/e11a;
e11b = E3(1,1);
e12b = E3(1,2);
e13b = E3(1,3);
d22b = S220-E3(2,2);
d23b = S121-E3(2,3);
d33b = S022-E3(3,3);
S220b = d22b+e13b^2/e11b;
S022b = d33b+e12b^2/e11b;
e11c = E6(1,1);
e12c = E6(1,2);
e13c = E6(1,3);
d22c = S022-E6(2,2);
d23c = S112-E6(2,3);
d33c = S202-E6(3,3);
S022c = d22c+e13c^2/e11c;
S202c = d33c+e12c^2/e11c;
% find minimum values
S220 = max([S220 S220a S220b]);
S202 = max([S202 S202a S202c]);
S022 = max([S022 S022b S022c]);
% recompute with updated values
e22a = S220-d22a;
e33a = S202-d33a;
e22b = S220-d22b;
e33b = S022-d33b;
e22c = S022-d22c;
e33c = S202-d33c;
%
S211 = realizability_S211(e11a,e22a,e33a,e12a,e13a,d23a,S211);
S121 = realizability_S211(e11b,e22b,e33b,e12b,e13b,d23b,S121);
S112 = realizability_S211(e11c,e22c,e33c,e12c,e13c,d23c,S112);
%




s220min = max([S220 S110^2 1-sqrt((S400-1)*(S040-1)) d22+(S301-d13)^2/(S400-d11) d22+(S031-d45)^2/(S040-d44)]);
s202min = max([S202 S101^2 1-sqrt((S400-1)*(S004-1)) d33+(S310-d12)^2/(S400-d11) d33+(S013-d56)^2/(S004-d66)]);
s022min = max([S022 S011^2 1-sqrt((S040-1)*(S004-1)) d55+(S103-d36)^2/(S004-d66) d55+(S130-d24)^2/(S040-d44)]);

X0 = S220-d22;
Y0 = S202-d33;
Z0 = S022-d55;


Xmax = s220max-d22;
Ymax = s202max-d33;
Zmax = s022max-d55;
Xmin = 0;
Ymin = 0;
Zmin = 0;

dDel1 = -(S011^2 - 2*S011*S101*S110 + S101^2 + S110^2 - 1)

[E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
% [E0a,E0b,E0c,E0d,E0e,E0f] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,s220,S030,S130,S040,...
%                       S101,S201,S301,S102,s202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,s022);
%
%dE4 = ([det(E1(1:4,1:4)) det(E2(1:4,1:4)) det(E3(1:4,1:4)) det(E4(1:4,1:4)) det(E5(1:4,1:4)) det(E6(1:4,1:4))])
%dE5 = ([det(E1(1:5,1:5)) det(E2(1:5,1:5)) det(E3(1:5,1:5)) det(E4(1:5,1:5)) det(E5(1:5,1:5)) det(E6(1:5,1:5))])
eE6 = eig(E1)
dE6 = det(E1)

% plots for 6 constraints
N = 200;
first_permutation
second_permutation
third_permutation
fourth_permutation
fifth_permutation
sixth_permutation
determinate6

% small = 1.d-6;
% xlim([small Xmax-small])
% ylim([small Ymax-small])
% zlim([small Zmax-small])

plot3(X0,Y0,Z0,'ko','MarkerFaceColor','k')

xlim([0.5*X0 1.2*X0])
ylim([0.1*Y0 3.*Y0])
zlim([0.1*Z0 3.*Z0])


