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

% Flux moment subscripts (i,j,k) for extracting from M5 array
% Fx moments: increment x index by 1
idx.Fx_subs = [
    2 1 1; 3 1 1; 4 1 1; 5 1 1; 6 1 1;  % M100-M500
    2 2 1; 3 2 1; 4 2 1; 5 2 1;        % M110-M410
    2 3 1; 3 3 1; 4 3 1;               % M120-M320
    2 4 1; 3 4 1; 2 5 1;               % M130,M230,M140
    2 1 2; 3 1 2; 4 1 2; 5 1 2;        % M101-M401
    2 1 3; 3 1 3; 4 1 3;               % M102-M302
    2 1 4; 3 1 4; 2 1 5;               % M103,M203,M104
    2 2 2; 3 2 2; 4 2 2;               % M111-M311
    2 3 2; 3 3 2; 2 4 2;               % M121,M221,M131
    2 2 3; 3 2 3; 2 2 4; 2 3 3         % M112,M212,M113,M122
];

% Fy moments: increment y index by 1
idx.Fy_subs = [
    1 2 1; 2 2 1; 3 2 1; 4 2 1; 5 2 1;  % M010-M410
    1 3 1; 2 3 1; 3 3 1; 4 3 1;        % M020-M320
    1 4 1; 2 4 1; 3 4 1;               % M030-M230
    1 5 1; 2 5 1; 1 6 1;               % M040,M140,M050
    1 2 2; 2 2 2; 3 2 2; 4 2 2;        % M011-M311
    1 2 3; 2 2 3; 3 2 3;               % M012-M212
    1 2 4; 2 2 4; 1 2 5;               % M013,M113,M014
    1 3 2; 2 3 2; 3 3 2;               % M021-M221
    1 4 2; 2 4 2; 1 5 2;               % M031,M131,M041
    1 3 3; 2 3 3; 1 3 4; 1 4 3         % M022,M122,M023,M032
];

% Fz moments: increment z index by 1
idx.Fz_subs = [
    1 1 2; 2 1 2; 3 1 2; 4 1 2; 5 1 2;  % M001-M401
    1 2 2; 2 2 2; 3 2 2; 4 2 2;        % M011-M311
    1 3 2; 2 3 2; 3 3 2;               % M021-M221
    1 4 2; 2 4 2; 1 5 2;               % M031,M131,M041
    1 1 3; 2 1 3; 3 1 3; 4 1 3;        % M002-M302
    1 1 4; 2 1 4; 3 1 4;               % M003-M203
    1 1 5; 2 1 5; 1 1 6;               % M004,M104,M005
    1 2 3; 2 2 3; 3 2 3;               % M012-M212
    1 3 3; 2 3 3; 1 4 3;               % M022,M122,M032
    1 2 4; 2 2 4; 1 2 5; 1 3 4         % M013,M113,M014,M023
];

end
