function M4 = pack_moments35(moments)
% PACK_MOMENTS35 Packs a moment structure into a 35-element vector
%
%   M4 = pack_moments35(moments) takes a structure with named moment fields
%   and returns a 35-element vector in the standard order
%
%   Input:
%       moments - Structure with 35 named moment fields
%   Output:
%       M4 - 35-element vector of moments

% Vector ordering:
% M = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
%      M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
%      M031,M012,M112,M013,M022]

M4 = zeros(35, 1);

% Order matters - this defines the canonical ordering
M4(1) = moments.M000;
M4(2) = moments.M100;
M4(3) = moments.M200;
M4(4) = moments.M300;
M4(5) = moments.M400;
M4(6) = moments.M010;
M4(7) = moments.M110;
M4(8) = moments.M210;
M4(9) = moments.M310;
M4(10) = moments.M020;
M4(11) = moments.M120;
M4(12) = moments.M220;
M4(13) = moments.M030;
M4(14) = moments.M130;
M4(15) = moments.M040;
M4(16) = moments.M001;
M4(17) = moments.M101;
M4(18) = moments.M201;
M4(19) = moments.M301;
M4(20) = moments.M002;
M4(21) = moments.M102;
M4(22) = moments.M202;
M4(23) = moments.M003;
M4(24) = moments.M103;
M4(25) = moments.M004;
M4(26) = moments.M011;
M4(27) = moments.M111;
M4(28) = moments.M211;
M4(29) = moments.M021;
M4(30) = moments.M121;
M4(31) = moments.M031;
M4(32) = moments.M012;
M4(33) = moments.M112;
M4(34) = moments.M013;
M4(35) = moments.M022;

end

