function [M5,C5,S5] = Moments5_3D(M4)
% Moments5_3D computes 3-D HyQMOM closure for Fluxes
%
%   Input:
% M4 = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%       M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%       M031,M012,M112,M013,M022]
% 
%   Output:
% M5 = [M500 M410 M320 M230 M140 M401 M302 M203 M104 M311 M221 M131 M212 M113 M122 M050 M041 M032 M023 M014 M005]
% C5 = [C500 C410 C320 C230 C140 C401 C302 C203 C104 C311 C221 C131 C212 C113 C122 C050 C041 C032 C023 C014 C005]
% S5 = [S500 S410 S320 S230 S140 S401 S302 S203 S104 S311 S221 S131 S212 S113 S122 S050 S041 S032 S023 S014 S005]

% Get mean velocities using helper
[umean, vmean, wmean] = means_from_M(M4);

% Compute central and standardized moments
[C4,S4] = M2CS4_35(M4);

% Get square roots of second-order central moments
[sC200, sC020, sC002] = sC_from_C4(C4);

% Call hyqmom_3D with standardized moments (reuse existing function)
[S500,S410,S320,S230,S140,S401,S302,S203,S104,S311,S221,S131,S212,S113,S122,S050,S041,S032,S023,S014,S005] = ...
    hyqmom_3D(S4(4),S4(5),S4(7),S4(8),S4(9),S4(11),S4(12),S4(13),S4(14),S4(15),...
              S4(17),S4(18),S4(19),S4(21),S4(22),S4(23),S4(24),S4(25),S4(26),...
              S4(27),S4(28),S4(29),S4(30),S4(31),S4(32),S4(33),S4(34),S4(35));

% Package S5 as cell array for scaling
S5_cell = {S500,S410,S320,S230,S140,S401,S302,S203,S104,S311,S221,S131,S212,S113,S122,S050,S041,S032,S023,S014,S005};

% Scale standardized moments to central moments using helper
C5_cell = scale_S5_to_C5(S5_cell, sC200, sC020, sC002);

% Convert to arrays for output
C5 = cell2mat(C5_cell);

% 5th-order moments from central moments
idx = moment_indices();
M000 = M4(idx.m000);
M5 = C5toM5_3D(M000,umean,vmean,wmean,C4(3),C4(7),C4(17),C4(10),C4(26),C4(20),...
               C4(4),C4(8),C4(18),C4(11),C4(27),C4(21),C4(13),C4(29),C4(32),C4(23),...
               C4(5),C4(9),C4(19),C4(12),C4(28),C4(22),C4(14),C4(30),C4(33),C4(24),...
               C4(15),C4(31),C4(35),C4(34),C4(25),...
               C5_cell{1},C5_cell{2},C5_cell{3},C5_cell{4},C5_cell{5},C5_cell{6},...
               C5_cell{7},C5_cell{8},C5_cell{9},C5_cell{10},C5_cell{11},C5_cell{12},...
               C5_cell{13},C5_cell{14},C5_cell{15},C5_cell{16},C5_cell{17},C5_cell{18},...
               C5_cell{19},C5_cell{20},C5_cell{21});
% Extract specific moments from M5 array for output
M500 = M5(6,1,1); M410 = M5(5,2,1); M320 = M5(4,3,1); M230 = M5(3,4,1); M140 = M5(2,5,1);
M401 = M5(5,1,2); M302 = M5(4,1,3); M203 = M5(3,1,4); M104 = M5(2,1,5);
M311 = M5(4,2,2); M221 = M5(3,3,2); M131 = M5(2,4,2);
M212 = M5(3,2,3); M113 = M5(2,2,4); M122 = M5(2,3,3);
M050 = M5(1,6,1); M041 = M5(1,5,2); M032 = M5(1,4,3); M023 = M5(1,3,4); M014 = M5(1,2,5); M005 = M5(1,1,6);

% Create output arrays
M5 = [M500 M410 M320 M230 M140 M401 M302 M203 M104 M311 M221 M131 M212 M113 M122 M050 M041 M032 M023 M014 M005]';
C5 = C5';  % Already in correct order from scale_S5_to_C5
S5 = cell2mat(S5_cell)';

end