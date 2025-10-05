function [C4,S4] = M2CS4_35(M4)
% M2CS4_35 Computes central and standardized moments from raw moments
%   [C4,S4] = M2CS4_35(M4) converts 35 raw moments to central (C4) and
%   standardized (S4) moments
% Input:  M4 - 35-element vector of raw moments
% Output: C4 - 35-element vector of central moments
%         S4 - 35-element vector of standardized moments

% Precompute linear indices (persistent = computed only once)
persistent idx_c;
if isempty(idx_c)
    % Manually computed sub2ind([5 5 5], i, j, k) = i + 5*(j-1) + 25*(k-1)
    idx_c = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 16, 17, 21, ...
             26, 27, 28, 29, 51, 52, 53, 76, 77, 101, 31, 32, 33, 36, 37, ...
             41, 56, 57, 81, 61];
end

C4 = zeros(size(M4));
S4 = zeros(size(M4));

% Extract moments (unpacking for readability in symbolic function call)
m = moment_struct('from_vector', M4);

% Compute central moments using symbolic function
C = M4toC4_3D(m.M000,m.M100,m.M200,m.M300,m.M400,m.M010,m.M110,m.M210,m.M310,...
              m.M020,m.M120,m.M220,m.M030,m.M130,m.M040,m.M001,m.M101,m.M201,...
              m.M301,m.M002,m.M102,m.M202,m.M003,m.M103,m.M004,m.M011,m.M111,...
              m.M211,m.M021,m.M121,m.M031,m.M012,m.M112,m.M013,m.M022);

% Extract variances for standardization (indices from 3D array)
sC200 = sqrt(max(C(3,1,1), eps));
sC020 = sqrt(max(C(1,3,1), eps));
sC002 = sqrt(max(C(1,1,3), eps));

% Compute standardized moments efficiently using array operations
% Zero-th to second order moments (fixed values)
S = zeros(5,5,5);
S(1,1,1) = 1;  % S000 = 1
S(3,1,1) = 1;  % S200 = 1
S(1,3,1) = 1;  % S020 = 1
S(1,1,3) = 1;  % S002 = 1

% Cross-correlations (second order)
S(2,2,1) = C(2,2,1)/(sC200*sC020);  % S110
S(2,1,2) = C(2,1,2)/(sC200*sC002);  % S101
S(1,2,2) = C(1,2,2)/(sC020*sC002);  % S011

% Third order moments
S(4,1,1) = C(4,1,1)/sC200^3;  % S300
S(3,2,1) = C(3,2,1)/(sC200^2*sC020);  % S210
S(3,1,2) = C(3,1,2)/(sC200^2*sC002);  % S201
S(2,3,1) = C(2,3,1)/(sC200*sC020^2);  % S120
S(2,2,2) = C(2,2,2)/(sC200*sC020*sC002);  % S111
S(2,1,3) = C(2,1,3)/(sC200*sC002^2);  % S102
S(1,4,1) = C(1,4,1)/sC020^3;  % S030
S(1,3,2) = C(1,3,2)/(sC020^2*sC002);  % S021
S(1,2,3) = C(1,2,3)/(sC020*sC002^2);  % S012
S(1,1,4) = C(1,1,4)/sC002^3;  % S003

% Fourth order moments
S(5,1,1) = C(5,1,1)/sC200^4;  % S400
S(4,2,1) = C(4,2,1)/(sC200^3*sC020);  % S310
S(4,1,2) = C(4,1,2)/(sC200^3*sC002);  % S301
S(3,3,1) = C(3,3,1)/(sC200^2*sC020^2);  % S220
S(3,2,2) = C(3,2,2)/(sC200^2*sC020*sC002);  % S211
S(3,1,3) = C(3,1,3)/(sC200^2*sC002^2);  % S202
S(2,4,1) = C(2,4,1)/(sC200*sC020^3);  % S130
S(2,3,2) = C(2,3,2)/(sC200*sC020^2*sC002);  % S121
S(2,2,3) = C(2,2,3)/(sC200*sC020*sC002^2);  % S112
S(2,1,4) = C(2,1,4)/(sC200*sC002^3);  % S103
S(1,5,1) = C(1,5,1)/sC020^4;  % S040
S(1,4,2) = C(1,4,2)/(sC020^3*sC002);  % S031
S(1,3,3) = C(1,3,3)/(sC020^2*sC002^2);  % S022
S(1,2,4) = C(1,2,4)/(sC020*sC002^3);  % S013
S(1,1,5) = C(1,1,5)/sC002^4;  % S004

C4 = C(idx_c);
S4 = S(idx_c);

end