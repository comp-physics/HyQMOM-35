% realizability regions for S220, S202, S002
clear
clc
close all

% select univariate moments
S300=0.2;
S030=S300;
S003=-0.1;
S400=4;
S040=3;
S004=3;
%
H200=S400-S300^2-1
H020=S040-S030^2-1
H002=S004-S003^2-1
% initial guesses based on absolute max
s220max = 1 + sqrt((S400-1)*(S040-1));
s202max = 1 + sqrt((S400-1)*(S004-1));
s022max = 1 + sqrt((S040-1)*(S004-1));
%
S220 = 0.95*s220max
S202 = 0.2*s202max
S022 = 0.2*s022max
%
% select covariances
S110= 0.9999;
S101= 0.;
S011= 0.;

% select third-order moments
S210=0*S110*S300;
S201=0;
S120=0*S110*S030;
S102=0;%S101*S003;
S021=0;%S011*S030;
S012=0;%S011*S003;
S111 = 0;
% select fourth-order moments
S310=0.;%S110*S400;
S301=0.;%S101*S400;
S130=0.;%S110*S040;
S103=0.;%S101*S004;
S031=0.;%S011*S040;
S013=0.;%S011*S004;
%
S211=0;
S121=0;
S112=0;
%
% % check realizability of selected moments
% [S110,S210,S310,S120,S220,S130,S101,S201,S301,S102,S202,S103,...
%           S011,S111,S211,S021,S121,S031,S012,S112,S013,S022] = ...
% realizable_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
%                  S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
%                  S211,S021,S121,S031,S012,S112,S013,S022);
% check realizability of selected moments
[S110,S210,S310,S120,S220,S130,S101,S201,S301,S102,S202,S103,...
          S011,S111,S211,S021,S121,S031,S012,S112,S013,S022] = ...
realizable_3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                 S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,...
                 S211,S021,S121,S031,S012,S112,S013,S022);
% find D matrix
%
Del1 = [1 0    0    0;...
        0 1    S110 S101;...
        0 S110 1    S011;...
        0 S101 S011 1];
%
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
%
s220min = max([0 S110^2 1-sqrt((S400-1)*(S040-1)) ]);
s202min = max([0 S101^2 1-sqrt((S400-1)*(S004-1)) ]);
s022min = max([0 S011^2 1-sqrt((S040-1)*(S004-1)) ]);
%
X0 = S220-d22
Y0 = S202-d33
Z0 = S022-d55
%
Xmax = s220max-d22;
Ymax = s202max-d33;
Zmax = s022max-d55;
Xmin = s220min-d22;
Ymin = s202min-d33;
Zmin = s022min-d55;

dDel1 = -(S011^2 - 2*S011*S101*S110 + S101^2 + S110^2 - 1)

[E1,E2,E3,E4,E5,E6] = delta2star3D_permutation(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,...
                      S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022);
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

small = 1.d-3;
xlim([Xmin+small Xmax-small])
ylim([Ymin+small Ymax-small])
zlim([Zmin+small Zmax-small])

plot3(X0,Y0,Z0,'ko','MarkerFaceColor','k')

xlim([0.5*X0 1.1*X0])
ylim([0.1*Y0 3.*Y0])
zlim([0.1*Z0 3.*Z0])


