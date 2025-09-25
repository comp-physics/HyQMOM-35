function idx = moment_indices()
% moment_indices Create index map for 35 3D moments to avoid magic numbers
%
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]

idx = struct();

% Commonly used moment indices
idx.m000 = 1;   % density
idx.m100 = 2;   % x-momentum
idx.m010 = 6;   % y-momentum  
idx.m001 = 16;  % z-momentum

% Second-order moments
idx.m200 = 3;
idx.m020 = 10;
idx.m002 = 20;
idx.m110 = 7;
idx.m101 = 17;
idx.m011 = 26;

% Third-order moments
idx.m300 = 4;
idx.m030 = 13;
idx.m003 = 23;

% Fourth-order moments  
idx.m400 = 5;
idx.m040 = 15;
idx.m004 = 25;

% Central moment indices (for C array)
idx.C200 = 3;
idx.C020 = 10; 
idx.C002 = 20;

% Standardized moment indices (for S array)
idx.S300 = 4;
idx.S030 = 13;
idx.S003 = 23;
idx.S400 = 5;
idx.S040 = 15;
idx.S004 = 25;

% 1D moment sequences for eigenvalue computation
idx.x_moments = [1, 2, 3, 4, 5];      % m000, m100, m200, m300, m400
idx.y_moments = [1, 6, 10, 13, 15];   % m000, m010, m020, m030, m040

end
