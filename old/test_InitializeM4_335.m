% test function M = InitializeM4_35(M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002) 
clear 
clc

% initial densities
rho = 1;
% initial mean velocities
U0 = 0;
V0 = 1;
W0 = 0;
% initial correlation coefficients for joint Gaussian
r110 = -1;
r101 = -1;
r011 = 1;
% set initial conditions to joint Gaussian with covariance
C200 = 1;
C020 = 1;
C002 = 1;
C110 = r110*sqrt(C200*C020);
C101 = r101*sqrt(C200*C002);
C011 = r011*sqrt(C020*C002);
% initialize moments
M4 = InitializeM4_35(rho,U0,V0,W0,C200,C110,C101,C020,C011,C002)
%
[C4,S4] = M2CS4_35(M4)
[MM5,CC5,SS5] = Moments5_3D(M4)
