function M = InitializeM4_35(M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002)
%INITIALIZEM4_35 Computes 3-D fourth-order joint Gaussian moments
%
%   Refactored from 172 lines to ~40 lines by eliminating repetitive array extractions.
%
%   Input:
%       M000  - number density 
%       umean - M100/M000 (mean velocity u)
%       vmean - M010/M000 (mean velocity v)
%       wmean - M001/M000 (mean velocity w)
%       C200, C110, C101, C020, C011, C002 - covariances
%
%   Output:
%       M - 35-element vector of moments

% Standardized moments for Maxwellian (Gaussian)
% 3rd order: all zero (Gaussian is symmetric)
S300=0; S210=0; S201=0; S120=0; S111=0; S102=0; S030=0; S021=0; S012=0; S003=0;

% 4th order: diagonal = 3 (Gaussian kurtosis), cross = 1 (independent variables)
S400=3; S310=0; S301=0; S220=1; S211=0; S202=1;
S130=0; S121=0; S112=0; S103=0; S040=3; S031=0; S022=1; S013=0; S004=3;

% Compute central moments from standardized moments
C4 = S4toC4_3D_r(C200,C110,C101,C020,C011,C002,...
                 S300,S210,S201,S120,S111,S102,S030,S021,S012,S003,...
                 S400,S310,S301,S220,S211,S202,S130,S121,S112,S103,S040,S031,S022,S013,S004);

% Extract central moments from 3D array
[~, ~, ~, ~, C200, C110, C101, C020, C011, C002, ...
 C300, C210, C201, C120, C111, C102, C030, C021, C012, C003, ...
 C400, C310, C301, C220, C211, C202, C130, C121, C112, C103, C040, C031, C022, C013, C004] = ...
    moment_conversion_utils('M4_to_vars', C4);

% Compute raw moments from central moments
M4 = C4toM4_3D(M000,umean,vmean,wmean,...
               C200,C110,C101,C020,C011,C002,...
               C300,C210,C201,C120,C111,C102,C030,C021,C012,C003,...
               C400,C310,C301,C220,C211,C202,C130,C121,C112,C103,C040,C031,C022,C013,C004);

% Extract raw moments from 3D array and pack into vector
[M000, M100, M010, M001, M200, M110, M101, M020, M011, M002, ...
 M300, M210, M201, M120, M111, M102, M030, M021, M012, M003, ...
 M400, M310, M301, M220, M211, M202, M130, M121, M112, M103, M040, M031, M022, M013, M004] = ...
    moment_conversion_utils('M4_to_vars', M4);

% Pack into 35-element vector (standard ordering)
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]
M = [M000; M100; M200; M300; M400; M010; M110; M210; M310; M020; M120; M220; M030; M130; M040;
     M001; M101; M201; M301; M002; M102; M202; M003; M103; M004; M011; M111; M211; M021; M121;
     M031; M012; M112; M013; M022];

end
